//Copyright 2017 Laurent Monasse

/*
  This file is part of Mka3D.
  
  Mka3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Mka3D is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Mka3D.  If not, see <http://www.gnu.org/licenses/>.
*/

/*!
\authors Laurent Monasse and Adela Puscas
 *  \file solide.cpp
 *  \brief Definition of the methods for solid classes.
 * Procedures specific to coupling are preceded by a "warning".
 */
#include "solide.hpp"
#include "particule.hpp"
#include "face.hpp"
#include "vertex.hpp"
#include "geometry.hpp"
#include "vitesse.hpp"
#include "forces_ext.hpp"
#include <iostream>
#include <map>
#include <set>
#include <string>
#ifndef SOLIDE_CPP
#define SOLIDE_CPP

Solide::Solide(const double& E, const double& nu){
  lambda = E * nu / (1.+nu) / (1. - 2.*nu);
  mu = E / 2. / (1.+nu);
}

Solide::Solide(){
  lambda = 0.;
  mu = 0.;
}

Solide::~Solide(){   
}

Solide & Solide::operator=(const Solide &S){
  assert(this != &S);
  //solide.resize(S.solide.size());
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    solide[P->first]= P->second;
  }
}

void Solide::Init(const char* s1, const char* s2, const char* s3, const char* s4, const bool& rep, const int& numrep, const double& rho){
  std::ifstream noeuds(s1,ios::in);
  std::ifstream elements(s2,ios::in);
  std::ifstream voisins(s3,ios::in);
  std::ifstream import_faces(s4,ios::in);
  if(not(noeuds) || not(elements) || not(voisins) || not(import_faces))
    cout <<"ouverture du maillage ratee" << endl;

  //Importation des vertex
  string s;
  char aux;
  string ligne;
  int nbr_vertex;
  getline(noeuds, ligne); //Lecture de la première ligne du fichier des vertex
  istringstream  stdd(ligne);
  stdd >> nbr_vertex; //Nombre de vertex

  for(int i=0; i < nbr_vertex ; i++) {
    getline(noeuds, ligne);
    istringstream  stm(ligne);
    int id; //Numéro du vertex + 1
    double x,y,z;
    stm >> id >> x >> y >> z;
    vertex.push_back(Vertex(Point_3(x,y,z))); //Vertex sont donnés dans l'ordre
  }

  //Importation des Particules
  getline(elements, ligne);
  istringstream stde(ligne);
  int nbr_elements;
  stde >> nbr_elements; //On stocke le nombre d'éléments
  for(int i=0; i < nbr_elements ; i++) {
    getline(elements, ligne);
    istringstream  stm(ligne);
    int id;
    stm >> id;
    int v1,v2,v3,v4;
    stm >> v1 >> v2 >> v3 >> v4;
    //Ajout des vertex de la particule
    Particule p(id);
    p.vertices.push_back(v1 - 1);
    p.vertices.push_back(v2 - 1);
    p.vertices.push_back(v3 - 1);
    p.vertices.push_back(v4 - 1);

    //Calcul des quantités volumiques (liées particule)
    p.barycentre(); //Calcul du barycentre
    p.volume(); //calcul du volume
    p.m = rho * p.V;

    //Ajout de la particule dans le solide
    solide[id] = p;
  }
  
  //Comment fixer BC ? Voir comment les gérer en tetgen ! Fichier en plus surement

  //Importation des faces
  getline(import_faces, ligne);
  istringstream stdf(ligne);
  int nbr_faces;
  stdf >> nbr_faces;
  for(int i=0; i<nbr_faces ; i++) {
    getline(elements, ligne);
    istringstream  stm(ligne);
    int id;
    int v1,v2,v3;
    stm >> id >> v1 >> v2 >> v3; //Numéro de la face et des vertex + 1
    Face f();
    f.id = id - 1;
    f.nb_vertex = 3;
    f.vertex.push_back(v1 - 1); //Ajout du numéro des vertex
    f.vertex.push_back(v2 - 1);
    f.vertex.push_back(v3 - 1);
    faces.push_back(f);
  }

  //Importation des voisins
  getline(voisins, ligne);
  istringstream stdv(ligne);
  int nbr_elements_bis;
  stdv >> nbr_elements_bis; //On stocke le nombre d'éléments
  if(nbr_elements != nb_elements_bis)
    cout << "Problème avec la connectivité du maillage !" << endl;
  for(int i=0 ; i<nbr_elements ; i++) {
    getline(voisins, ligne);
    istringstream  stm(ligne);
    int id; //Numéro de la particule dont on récupère les voisins
    stm >> id;
    int v1,v2,v3,v4;
    stm >> v1 >> v2 >> v3 >> v4; //Id des voisins dans solide. -1 indique une particule sans voisins donc sur le bord

    //Ajout des voisins
    //face 1
    if(v1 >= 0) { //Sinon particule sur le bord et pas de voisin
      bool existence = false; //Test si face et voisins déjà liés
      for(int voi=0; voi<solide[v1].faces.size ; voi++) { //Regarde si dans les faces liées à la particule voisine, la particule apparaît comme voisine
	if(((solide[v1].faces[voi])->voisins)[1] == id)
	  existence = true; //Dans ce cas face déjà liée et voisins ok
	else if(((solide[v1].faces[voi])->voisins)[0] == id)
	  existence = true; //idem
      }

      if(not(existence)) {
	int id_face = 0; //Trouver le numéro de la face
	std::vector<int> vert_face(vertex_face(id, v1)); //Trouve quels sont les vertex de la face
	for(int j=0; j <faces.size ; j++) { //Trouve la face dans l'ensemble des faces
	  Face f();
	  f.vertex.push_back(solide[id].vertices[vert_face[0]]);
	  f.vertex.push_back(solide[id].vertices[vert_face[1]]);
	  f.vertex.push_back(solide[id].vertices[vert_face[2]]);
	  if(f == faces[j])
	    id_face = j;
	}

	solide[id].faces.push_back(id_face); //On ajoute la face dans la liste pour la particule
	solide[v1].faces.push_back(id_face); //idem pour la particule voisine
	faces[id_face].voisins.push_back(id); //Ajout de la particule dans la liste des voisins de la particule
	faces[id_face].voisins.push_back(v1); //Ajout de la particule voisine dans la liste des voisins de la particule     
	faces[id_face].comp_quantities(solide[id].vertices[vert_face[0]].pos, solide[id].vertices[vert_face[1]].pos, solide[id].vertices[vert_face[0]].pos, solide[id].vertices[vert_face[3]].pos); //Calcul de la normale sortante, surface et barycentre face
      }
    }

    //face 2
    if(v2 >= 0) { //Sinon particule sur le bord et pas de voisin
      existence = false; //Test si face et voisins déjà liés
      for(int voi=0; voi<solide[v2].faces.size ; voi++) { //Regarde si dans les faces liées à la particule voisine, la particule apparaît comme voisine
	if(((solide[v2].faces[voi])->voisins)[1] == id)
	  existence = true; //Dans ce cas face déjà liée et voisins ok
	else if(((solide[v2].faces[voi])->voisins)[0] == id)
	  existence = true; //idem
      }

      if(not(existence)) {
	int id_face = 0; //Trouver le numéro de la face
	std::vector<int> vert_face(vertex_face(id, v2)); //Trouve quels sont les vertex de la face
	for(int j=0; j <faces.size ; j++) { //Trouve la face dans l'ensemble des faces
	  Face f();
	  f.vertex.push_back(solide[id].vertices[vert_face[0]]);
	  f.vertex.push_back(solide[id].vertices[vert_face[1]]);
	  f.vertex.push_back(solide[id].vertices[vert_face[2]]);
	  if(f == faces[j])
	    id_face = j;
	}

	solide[id].faces.push_back(id_face); //On ajoute la face dans la liste pour la particule
	solide[v2].faces.push_back(id_face); //idem pour la particule voisine
	faces[id_face].voisins.push_back(id); //Ajout de la particule dans la liste des voisins de la particule
	faces[id_face].voisins.push_back(v2); //Ajout de la particule voisine dans la liste des voisins de la particule     
	faces[id_face].comp_quantities(solide[id].vertices[vert_face[0]].pos, solide[id].vertices[vert_face[1]].pos, solide[id].vertices[vert_face[0]].pos, solide[id].vertices[vert_face[3]].pos); //Calcul de la normale sortante, surface et barycentre face
      }
    }

    //face 3
    if(v3 >= 0) { //Sinon particule sur le bord et pas de voisin
      existence = false; //Test si face et voisins déjà liés
      for(int voi=0; voi<solide[v3].faces.size ; voi++) { //Regarde si dans les faces liées à la particule voisine, la particule apparaît comme voisine
	if(((solide[v3].faces[voi])->voisins)[1] == id)
	  existence = true; //Dans ce cas face déjà liée et voisins ok
	else if(((solide[v3].faces[voi])->voisins)[0] == id)
	  existence = true; //idem
      }

      if(not(existence)) {
	int id_face = 0; //Trouver le numéro de la face
	std::vector<int> vert_face(vertex_face(id, v3)); //Trouve quels sont les vertex de la face
	for(int j=0; j <faces.size ; j++) { //Trouve la face dans l'ensemble des faces
	  Face f();
	  f.vertex.push_back(solide[id].vertices[vert_face[0]]);
	  f.vertex.push_back(solide[id].vertices[vert_face[1]]);
	  f.vertex.push_back(solide[id].vertices[vert_face[2]]);
	  if(f == faces[j])
	    id_face = j;
	}

	solide[id].faces.push_back(id_face); //On ajoute la face dans la liste pour la particule
	solide[v3].faces.push_back(id_face); //idem pour la particule voisine
	faces[id_face].voisins.push_back(id); //Ajout de la particule dans la liste des voisins de la particule
	faces[id_face].voisins.push_back(v3); //Ajout de la particule voisine dans la liste des voisins de la particule     
	faces[id_face].comp_quantities(solide[id].vertices[vert_face[0]].pos, solide[id].vertices[vert_face[1]].pos, solide[id].vertices[vert_face[0]].pos, solide[id].vertices[vert_face[3]].pos); //Calcul de la normale sortante, surface et barycentre face
      }
    }

    //face 4
    if(v4 >= 0) { //Sinon particule sur le bord et pas de voisin
      existence = false; //Test si face et voisins déjà liés
      for(int voi=0; voi<solide[v4].faces.size ; voi++) { //Regarde si dans les faces liées à la particule voisine, la particule apparaît comme voisine
	if(((solide[v4].faces[voi])->voisins)[1] == id)
	  existence = true; //Dans ce cas face déjà liée et voisins ok
	else if(((solide[v4].faces[voi])->voisins)[0] == id)
	  existence = true; //idem
      }

      if(not(existence)) {
	int id_face = 0; //Trouver le numéro de la face
	std::vector<int> vert_face(vertex_face(id, v1)); //Trouve quels sont les vertex de la face
	for(int j=0; j <faces.size ; j++) { //Trouve la face dans l'ensemble des faces
	  Face f();
	  f.vertex.push_back(solide[id].vertices[vert_face[0]]);
	  f.vertex.push_back(solide[id].vertices[vert_face[1]]);
	  f.vertex.push_back(solide[id].vertices[vert_face[2]]);
	  if(f == faces[j])
	    id_face = j;
	}

	solide[id].faces.push_back(id_face); //On ajoute la face dans la liste pour la particule
	solide[v4].faces.push_back(id_face); //idem pour la particule voisine
	faces[id_face].voisins.push_back(id); //Ajout de la particule dans la liste des voisins de la particule
	faces[id_face].voisins.push_back(v4); //Ajout de la particule voisine dans la liste des voisins de la particule     
	faces[id_face].comp_quantities(solide[id].vertices[vert_face[0]].pos, solide[id].vertices[vert_face[1]].pos, solide[id].vertices[vert_face[0]].pos, solide[id].vertices[vert_face[3]].pos); //Calcul de la normale sortante, surface et barycentre face
      }
    }
  } //fin boucle importation des connectivités
  
  //Calcul du tetrahèdre associé à chaque face pour le calcul du gradient
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      if((F->voisins).size() != 0) {//Face pas au bord (car I_Dx = vec null dans ce cas...)
	//Recherche des particules pour avoir le tetra !
	//for(std::vector<Face>::iterator G=solide[(F->voisins)[0]].faces.begin();G!=solide[(F->voisins)[0]].faces.end();G++) {
	Particule *part(&solide[(F->voisins)[0]]); //Pointeur sur la particule voisine

	//Vérifier si on bien le tétra et ajouter les 2 valeurs nécessaires.
	(F->voisins).push_back((part->faces)[0]
	
      }
    }
  }
  
  Vector_3 aux(solide[(it->parts)[0]], x0);
  double D0 = sqrt(aux.squared_length()); //Distance en config initiale entre particules voisines
  double D_face = abs(Vector_3(x0, v2.pos) * face1.normale); //Distance centre particule-face
  //Calculer la suite ! Mettre plus haut ?
  }
}

std::vector<int> Solide::vertex_face(const int& particule, const int& voisin) {
  std::vector<int> resultat;
  for(int j=0; j <solide[voisin].vertices.size ; j++) {
    if(solide[particule].vertices[0] == solide[voisin].vertices[j])
      resultat.push_back(0);
  }
  for(int j=0; j <solide[voisin].vertices.size ; j++) {
    if(solide[particule].vertices[1] == solide[voisin].vertices[j])
      resultat.push_back(1);
  }
  for(int j=0; j <solide[voisin].vertices.size ; j++) {
    if(solide[particule].vertices[2] == solide[voisin].vertices[j])
      resultat.push_back(2);
  }
  for(int j=0; j <solide[voisin].vertices.size ; j++) {
    if(solide[particule].vertices[4] == solide[voisin].vertices[j])
      resultat.push_back(4);
  }

  for(int j=0; j <solide[particule].vertices.size ; j++) {
    if(solide[particule].vertices[j] != resultat[0] && solide[particule].vertices[j] != resultat[1] && solide[particule].vertices[j] != resultat[2])
      resultat.push_back(j); //Dernière particule hors de la face mais sert pour le calcul de la normale
  }
  return resultat;
}

void Solide::Solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).solve_position(dt, flag_2d, t, T);
  }
}

void Solide::Solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).solve_vitesse(dt, flag_2d, Amort, t , T);
  }
}

void Solide::Forces(const int& N_dim, const double& dt, const double& t, const double& T){
  Forces_internes(N_dim,nu,E, dt);
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      (P->second).Fi = (P->second).Fi + Forces_externes(t,T, *F);
  }
}

void Solide::stresses(const double& dt){ //Calcul de la contrainte dans chaque particule
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    (P->second).discrete_gradient.col2 = Vector_3(0., 0., 0.);
    (P->second).discrete_gradient.col3 = Vector_3(0., 0., 0.);
    for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++){
      if(F->voisin>=0){
	int part = F->voisin;
	Vector_3 nIJ = (*F).normale;
	Matrix Dij_n(tens_sym(solide[part].Dx + solide[part].u * dt/2. - (P->second).Dx - (P->second).u * dt/2.,  nIJ) ); //Quadrature au point milieu pour calcul des forces !
	(P->second).discrete_gradient += (*F).S / 2. * Dij_n / (P->second).V;

      }
    }
    //cout << "Trace dev Def : " << (((P->second).discrete_gradient).dev()).tr() << endl;
    (P->second).contrainte = lambda * ((P->second).discrete_gradient - (P->second).epsilon_p).tr() * unit() + 2*mu * ((P->second).discrete_gradient - (P->second).epsilon_p);
    //cout << "Trace dev Contrainte : " << (((P->second).contrainte).dev()).tr() << endl;

    bool plastifie = false;
    //Mettre ces valeurs dans le param.dat !!!!!
    double B = 292000000.; //En Pa. JC.
    double n = .31; //JC.
    double A = 90000000.; //En Pa. Vient de JC
    double H = 0.;
	
    (P->second).seuil_elas = A; // + B * pow((P->second).def_plas_cumulee, n);

    if(((P->second).contrainte - H * (P->second).epsilon_p).VM() > (P->second).seuil_elas) { //On sort du domaine élastique.
      plastifie = true;
      //Matrix n_elas(((*P).contrainte).dev() / (((*P).contrainte).dev()).norme() ); //Normale au domaine élastique de Von Mises
      Matrix n_elas( 1. / (((P->second).contrainte).dev()).norme() * ((P->second).contrainte).dev() ); //Normale au domaine élastique de Von Mises
      /*if((*P).n_elas_prev == -n_elas)
	cout << "Chargement dans sens oppose !" << endl;
      (*P).n_elas_prev = n_elas;*/
      //cout << "Trace n_elas : " << n_elas.tr() << endl;
      //cout << "Norme n_elas : " << n_elas.norme() << endl;
      //double delta_p = pow(((*P).contrainte.VM() - A) / B, 1./n) - (*P).def_plas_cumulee;
      double delta_p = (((P->second).contrainte - H * (P->second).epsilon_p).VM() - A) / (2*mu + H);
      (*P).def_plas_cumulee += delta_p;
      //cout << "Def plastique cumulee : " << (*P).def_plas_cumulee << endl;
      (P->second).epsilon_p += delta_p * n_elas;
      //cout << "Trace def plas : " << ((P->second).epsilon_p).tr() << endl; //Pb ! Non-nulle !!!!
      //cout << "Norme def plas : " << ((P->second).epsilon_p).norme() << endl;
      
      //((P->second).contrainte.VM() - A) / (2*mu) * n_elas;  //* signe( (P->second).contrainte ); //Plasticité parfaite
      //(P->second).contrainte = A * signe( (P->second).contrainte );
    }
  }

  if(plastifie)
    cout << "Plastification dans ce pas de temps !" << endl;
}


void Solide::Forces_internes(const double& dt){ //Calcul des forces pour chaque particule
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).Fi = Vector_3(0.,0.,0.);
    for(std::vector<Face *>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++){
      if((*F).voisin>=0){
	int part = (*F).voisin;
	Vector_3 nIJ = (*F).normale;
        Vector_3 Fij_elas( (*F).S / 2. * ( ((P->second).contrainte + solide[part].contrainte).tr() / 2. * nIJ + ((P->second).contrainte + solide[part].contrainte) / 2. * nIJ ) ); //Force du lien IJ !
	//cout << "Force : " << Fij_elas << endl;

	(P->second).Fi = (P->second).Fi + Fij_elas; // * nIJ; //Force sur particule
	
      }
      //cout << "Force interne : " << (P->second).Fi << endl;
    }
  }
}

double Solide::Energie(const int& N_dim, const double& nu, const double& E){
  return Energie_cinetique()+Energie_potentielle(N_dim, nu, E);
}

double Solide::Energie_cinetique(){
  double E = 0.;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    double u2 = ((P->second).u.squared_length());
    E += 1./2. * (P->second).m * u2;
  }
  //cout << "Energie cinetique : " << E << endl;
  return E;
}

double Solide::Energie_potentielle(const int& N_dim, const double& nu, const double& E){
  double Ep = 0.;

  double B = 292000000.; //En Pa. JC.
  double n = .31; //JC.
  double A = 90000000.; //En Pa. Vient de JC

  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    Ep += ( 0.5 * contraction_double(((P->second)).contrainte, ((P->second)).discrete_gradient - ((P->second)).epsilon_p)  + B * pow(((P->second)).def_plas_cumulee, 1. + n) / (n + 1.) + A * ((P->second)).def_plas_cumulee ) * ((P->second)).V;
  }
  //cout << "Energie potentielle : " << Ep << endl;
  return Ep;
}

double Solide::pas_temps(const double& t, const double& T, const double& cfls, const double& E, const double& nu, const double& rhos){
  double eps = 1e-14;//std::numeric_limits<double>::epsilon();
  double dt = std::numeric_limits<double>::infinity();
  //Restriction CFL liee aux forces internes
  double cs = sqrt(E*(1.-nu)/rhos/(1.+nu)/(1.-2.*nu));
  //Calcul du rayon de la sphï¿½re inscrite
  double sigma = 100000.;
  for(int i=0;i<size();i++){
    for(int j=0;j<solide[i].faces.size();j++){
      sigma = min(sigma,solide[i].faces[j].D0);
    }
  }
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(int j=0;j<(P->second).faces.size();j++){
      if((P->second).faces[j].voisin>=0){
	dt = min(dt,cfls*(P->second).faces[j].D0/cs);
	/*double Imin = min(min((P->second).I[0],(P->second).I[1]),(P->second).I[2]);
	double S = (P->second).faces[j].S;
	double D0 = (P->second).faces[j].D0;
	double kappa = 1.;
	double alphan = (2.+2.*nu-kappa)*E/4./(1.+nu)/S*((P->second).faces[j].Is+(P->second).faces[j].It);
	double alphas = E/4./(1.+nu)/S*((2.+2.*nu+kappa)*(P->second).faces[j].Is-(2.+2.*nu-kappa)*(P->second).faces[j].It);
	double alphat = E/4./(1.+nu)/S*((2.+2.*nu+kappa)*(P->second).faces[j].It-(2.+2.*nu-kappa)*(P->second).faces[j].Is);
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphan));
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphas));
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphat));*/
      }
    }
  }
  dt = min(dt,T-t);
  return dt;
}

/*void Solide::update_triangles(){
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).triangles_prev = (P->second).triangles;
    (P->second).normales_prev = (P->second).normales;
    (P->second).fluide_prev = (P->second).fluide;
    for(int it=0;it<(P->second).triangles.size();it++){
      (P->second).Points_interface_prev[it] = (P->second).Points_interface[it];
      (P->second).Triangles_interface_prev[it] = (P->second).Triangles_interface[it];
      (P->second).Position_Triangles_interface_prev[it] = (P->second).Position_Triangles_interface[it];
      (P->second).Points_interface[it].erase((P->second).Points_interface[it].begin(),(P->second).Points_interface[it].end());
      (P->second).Triangles_interface[it].erase((P->second).Triangles_interface[it].begin(),(P->second).Triangles_interface[it].end());	(P->second).Position_Triangles_interface[it].erase((P->second).Position_Triangles_interface[it].begin(), (P->second).Position_Triangles_interface[it].end());
    }
    (P->second).triangles.erase((P->second).triangles.begin(),(P->second).triangles.end());
    (P->second).normales.erase((P->second).normales.begin(),(P->second).normales.end());
    (P->second).fluide.erase((P->second).fluide.begin(),(P->second).fluide.end());
    (P->second).vide.erase((P->second).vide.begin(),(P->second).vide.end());
		
    //Calcul de la nouvelle position des triangles
    for(int f=0;f<(P->second).faces.size();f++){
      Point_3 s,r,v,t;
			
      if((P->second).faces[f].size() == 3){
	vector<Point_3> ri,vi,si ;
	for(int part=0; part<(P->second).faces[f].vertex[0].size();part++){
	  int p = (P->second).faces[f].vertex[0].particules[part];
	  ri.push_back(solide[p].mvt_t((P->second).faces[f].vertex[0].pos));
	}
	r = centroid(ri.begin(),ri.end());
				
				
	for(int part=0;part<(P->second).faces[f].vertex[1].size();part++){
	  int p = (P->second).faces[f].vertex[1].particules[part];
	  vi.push_back(solide[p].mvt_t((P->second).faces[f].vertex[1].pos));
	}
	v = centroid(vi.begin(),vi.end());

	for(int part=0;part<(P->second).faces[f].vertex[2].size();part++){
	  int p = (P->second).faces[f].vertex[2].particules[part];
	  si.push_back(solide[p].mvt_t((P->second).faces[f].vertex[2].pos));
	}
	s = centroid(si.begin(),si.end());
				
	Vector_3 vect0(r,v);
	Vector_3 vect1(r,s);
	Triangle_3 Tri(r,v,s);
	(P->second).triangles.push_back(Tri);
	Vector_3 normale = cross_product(vect0,vect1);
	normale = normale*(1./sqrt((normale.squared_length())));
	(P->second).normales.push_back(normale);
	if((P->second).faces[f].voisin < 0){
	  (P->second).fluide.push_back(true);
	} else {
	  (P->second).fluide.push_back(false);
	}
	if( (P->second).faces[f].voisin == -2){
	  (P->second).vide.push_back(true);
	} else {
	  (P->second).vide.push_back(false);
	}
      }
      else{
			
	vector<Point_3> si;
	si.push_back((P->second).mvt_t((P->second).faces[f].centre));
	int j = (P->second).faces[f].voisin;
	if(j>=0){
	  si.push_back(solide[j].mvt_t((P->second).faces[f].centre));
	}
	s = centroid(si.begin(),si.end());
			
	for(int k=0;k<(P->second).faces[f].size();k++){
	  int kp = (k+1)%((P->second).faces[f].size());
	  vector<Point_3> ri,vi;
	  for(int part=0;part<(P->second).faces[f].vertex[k].size();part++){
	    int p = (P->second).faces[f].vertex[k].particules[part];
	    ri.push_back(solide[p].mvt_t((P->second).faces[f].vertex[k].pos));
	  }
	  r = centroid(ri.begin(),ri.end());
	  for(int part=0;part<(P->second).faces[f].vertex[kp].size();part++){
	    int p = (P->second).faces[f].vertex[kp].particules[part];
	    vi.push_back(solide[p].mvt_t((P->second).faces[f].vertex[kp].pos));
	  }
	  v = centroid(vi.begin(),vi.end());
	  Vector_3 vect0(s,r);
	  Vector_3 vect1(s,v);
	  Triangle_3 Tri(s,r,v);
	  (P->second).triangles.push_back(Tri);
	  Vector_3 normale = cross_product(vect0,vect1);
	  normale = normale*(1./sqrt((normale.squared_length())));			  
	  (P->second).normales.push_back(normale);
	  if((P->second).faces[f].voisin < 0){
	    (P->second).fluide.push_back(true);
	  } 
	  else {
	    (P->second).fluide.push_back(false);
	  }
	  if( (P->second).faces[f].voisin == -2){
	    (P->second).vide.push_back(true);
	  } else {
	    (P->second).vide.push_back(false);
	  }
	}
      }
			
    }//Calcul de la nouvelle position des triangles
  }
  }*/

void Solide::Impression(const int &n){ //Sortie au format vtk
  int nb_part = solide.size();
  std::ostringstream oss;
  oss << "solide" << n << ".vtk";
  string s = oss.str();
  const char* const solidevtk = s.c_str();
    
  //Ouverture des flux en donne en ecriture
  std::ofstream vtk;
  vtk.open(solidevtk,ios::out);
  if(not(vtk.is_open()))
    cout <<"ouverture de solide" << n << ".vtk rate" << endl;
  //vtk << setprecision(15);
  
  //Pour tetras !
  int nb_points = 4 * nb_part;
  int nb_faces = 4 * nb_part;
  int size = 4 * nb_faces;
  /*for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      size += F->nb_vertex; //Egale à nb_points au final ?
  }*/
  size += nb_faces;
    
  //Initialisation du fichier vtk
  vtk << "# vtk DataFile Version 3.0" << endl;
  vtk << "#Simulation Euler" << endl;
  vtk << "ASCII" << endl;
  vtk<< "\n";
  vtk << "DATASET UNSTRUCTURED_GRID" << endl;
  vtk << "POINTS " << nb_points << " DOUBLE" << endl;
    
  //Sortie des points
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Point_3>::iterator V=(P->second).vertices.begin();V!=(P->second).vertices.end();V++) {
      vtk << (P->second).mvt_t(*V) << endl;
    }
  }
  vtk << "\n";
    
  //Sortie des faces
  int point_tmp=0;
  vtk << "CELLS " << nb_faces << " " << size << endl;
  int compteur_vertex = 0;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face*>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      vtk << F->nb_vertex;
      for(int k=0 ; k<F->nb_vertex ; k++)
	vtk << " " << compteur_vertex + (F->vertex)[k];
      vtk << endl;
    }
    compteur_vertex += (P->second).vertices.size();
    //vtk << endl;
  }
  vtk << "\n";
  vtk << "CELL_TYPES " << nb_faces << endl;
  for(int i=0;i<nb_faces;i++){
    vtk << 7 << endl;
  }
  vtk << "\n";
  vtk << "CELL_DATA " << nb_faces << endl;
  //Deplacement
  vtk << "VECTORS deplacement double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face *>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      vtk << (P->second).Dx << endl;
  }
  vtk << "\n";
  //Vitesse
  vtk << "VECTORS vitesse double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face *>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      vtk << (P->second).u << endl;
  }
  vtk << "\n";
  //Contrainte
  vtk << "TENSORS contraintes double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face *>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      vtk << (P->second).contrainte.col1[0] << " " << (P->second).contrainte.col1[1] << " " << (P->second).contrainte.col1[2] << endl;
      vtk << (P->second).contrainte.col2[0] << " " << (P->second).contrainte.col2[1] << " " << (P->second).contrainte.col2[2] << endl;
      vtk << (P->second).contrainte.col3[0] << " " << (P->second).contrainte.col3[1] << " " << (P->second).contrainte.col3[2] << endl;
    }
  }
  vtk << "\n";
  //Déformations
  vtk << "TENSORS deformations double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face *>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      vtk << (P->second).discrete_gradient.col1[0] << " " << (P->second).discrete_gradient.col1[1] << " " << (P->second).discrete_gradient.col1[2] << endl;
      vtk << (P->second).discrete_gradient.col2[0] << " " << (P->second).discrete_gradient.col2[1] << " " << (P->second).discrete_gradient.col2[2] << endl;
      vtk << (P->second).discrete_gradient.col3[0] << " " << (P->second).discrete_gradient.col3[1] << " " << (P->second).discrete_gradient.col3[2] << endl;
    }
  }
  vtk << "\n";
  //Epsilon_p
  vtk << "TENSORS epsilon_p double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face *>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      vtk << (P->second).epsilon_p.col1[0] << " " << (P->second).epsilon_p.col1[1] << " " << (P->second).epsilon_p.col1[2] << endl;
      vtk << (P->second).epsilon_p.col2[0] << " " << (P->second).epsilon_p.col2[1] << " " << (P->second).epsilon_p.col2[2] << endl;
      vtk << (P->second).epsilon_p.col3[0] << " " << (P->second).epsilon_p.col3[1] << " " << (P->second).epsilon_p.col3[2] << endl;
    }
  }
  vtk << "\n";
  //Deformation plastique cumulée
  vtk << "SCALARS p double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face *>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      vtk << (P->second).def_plas_cumulee << endl;
  }
  vtk << "\n";
  vtk.close();
}

#endif
