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

Solide::Solide(const double& E, const double& nu, const double& B1, const double& n1, const double& A1, const double& H1){
  lambda = E * nu / (1.+nu) / (1. - 2.*nu);
  mu = E / 2. / (1.+nu);
  A = A1;
  B = B1;
  n = n1;
  H = H1;
}

Solide::Solide(){
  lambda = 0.;
  mu = 0.;
  A = 0.;
  B = 0.;
  n = 0.;
  H = 0.;
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

void Solide::Init(const char* s1, const char* s2, const char* s3, const bool& rep, const int& numrep, const double& rho){
  std::ifstream noeuds(s1,ios::in);
  std::ifstream elements(s2,ios::in);
  //std::ifstream voisins(s3,ios::in);
  std::ifstream import_faces(s3,ios::in);
  if(not(noeuds) || not(elements) /*|| not(voisins)*/ || not(import_faces))
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
    vertex.push_back(Vertex(Point_3(x,y,z), id - 1)); //Vertex sont donnés dans l'ordre
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
    p.vertices.push_back(vertex[v1 - 1].pos);
    p.vertices.push_back(vertex[v2 - 1].pos);
    p.vertices.push_back(vertex[v3 - 1].pos);
    p.vertices.push_back(vertex[v4 - 1].pos);

    //Calcul des quantités volumiques (liées particule)
    p.barycentre(); //Calcul du barycentre
    p.volume(); //calcul du volume
    p.m = rho * p.V;

    //Ajout de la particule dans le solide
    solide[id] = p;
  }
  
  //Comment fixer BC ? Voir comment les gérer en tetgen ! Fichier en plus surement

  //Importation des faces et des connectivités
  getline(import_faces, ligne);
  istringstream stdf(ligne);
  int nbr_faces;
  stdf >> nbr_faces;
  for(int i=0; i<nbr_faces ; i++) {
    getline(import_faces, ligne);
    istringstream  stm(ligne);
    int id;
    int v1,v2,v3; //Les vertex de la face
    int part_1, part_2; //Le numéro des particules voisines
    stm >> id >> v1 >> v2 >> v3 >> part_1 >> part_2; //Numéro de la face, des vertex + 1 et des voisins (bon numéro)
    Face f;
    f.id = id - 1;
    f.vertex.push_back(v1 - 1); //Ajout du numéro des vertex
    f.vertex.push_back(v2 - 1);
    f.vertex.push_back(v3 - 1);
    if(part_1 == -1 || part_2 == -1)
      f.BC = 1; //face au bord
    else {
      f.BC = 0; //Face pas au bord
      f.D0 = sqrt(Vector_3(solide[part_1].x0, solide[part_2].x0).squared_length());
    }
    if(part_1 >=0) { //Ajout du numéro des voisins dans la face
      f.voisins.push_back(part_1);
      f.voisins.push_back(part_2);
    }
    else { //La première valeur de voisin est la seule particule qui contient la face sur le bord
      f.voisins.push_back(part_2);
      f.voisins.push_back(part_1);
    }
    bool calcul_normales = false;
    if(part_1 >= 0) { //cad particule pas sur le bord
      solide[part_1].faces.push_back(f.id); //Ajout du numéro de la face dans la liste ds voisins de chaque particule
      f.comp_quantities(vertex[v1-1].pos, vertex[v2-1].pos, vertex[v3-1].pos); //Calcul de la normale sortante, surface et barycentre face
      calcul_normales = true;
      
    }
    if(part_2 >= 0) { //cad particule pas sur le bord
      solide[part_2].faces.push_back(f.id);
      if(not(calcul_normales)) {
	f.comp_quantities(vertex[v1-1].pos, vertex[v2-1].pos, vertex[v3-1].pos); //, solide[part_2].vertices[id_vertex_hors_face]); //Calcul de la normale sortante, surface et barycentre face
      }
    }
    //Vérification du sens de la normale
    if(part_1 != -1 && part_2 != -1) { //Face pas au bord
      if(Vector_3(solide[part_1].x0, solide[part_2].x0) * f.normale  < 0.)
	f.normale = -f.normale;
    }
    if(part_1 == -1 || part_2 == -1) { //Face au bord
      f.normale = -f.normale; //normale forcement dans mauvais sens
    }
    faces.push_back(f);
  }

  //Calcul du tetrahèdre associé à chaque face pour le calcul du gradient
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){ //Boucle sur toutes les faces
    int part_1 = F->voisins[0];
    int part_2 = F->voisins[1];
    bool tetra_ok = false;
    if(not(part_1 == -1)) { //Face dans le bulk et pas sur le bord
      int voisin1 = -1, voisin2 = -1; //Vont être les 2 autres particules composant le tetra associé à la face 
      //Boucle dans les faces de la première particule
      for(std::vector<int>::iterator G=solide[part_1].faces.begin();G!=solide[part_1].faces.end();G++) {
	if(not(faces[*G] == *F) && faces[*G].voisins[0] != part_1 && faces[*G].voisins[0] != -1 && voisin1 == -1)
	  voisin1 = faces[*G].voisins[0];
	else if(not(faces[*G] == *F) && faces[*G].voisins[1] != -1 && faces[*G].voisins[1] != part_1 && voisin1 == -1)
	  voisin1 = faces[*G].voisins[1];
      }
      for(std::vector<int>::iterator G=solide[part_1].faces.begin();G!=solide[part_1].faces.end();G++) {
	if(not(faces[*G] == *F) && faces[*G].voisins[0] != -1 && faces[*G].voisins[0] != part_1 && voisin1 != -1)
	  voisin2 = faces[*G].voisins[0];
	else if(not(faces[*G] == *F) && faces[*G].voisins[1] != -1 && faces[*G].voisins[1] != part_1 && voisin1 != -1)
	  voisin2 = faces[*G].voisins[1]; 

	//Vérifie que le centre de la face est dans le tetra et tetra pas aplati
	double vol = cross_product(Vector_3(solide[part_1].x0,solide[part_2].x0),Vector_3(solide[part_1].x0,solide[voisin1].x0))*Vector_3(solide[part_1].x0,solide[voisin2].x0)/6.; //Volume du tetra associé à la face
	double c_part_1 = cross_product(Vector_3(F->centre,solide[part_2].x0),Vector_3(F->centre,solide[voisin1].x0))*Vector_3(F->centre,solide[voisin2].x0)/6. / vol;
	double c_part_2 = cross_product(Vector_3(F->centre,solide[part_1].x0),Vector_3(F->centre,solide[voisin1].x0))*Vector_3(F->centre,solide[voisin2].x0)/6. / vol;
	double c_voisin1 = cross_product(Vector_3(F->centre,solide[part_1].x0),Vector_3(F->centre,solide[part_2].x0))*Vector_3(F->centre,solide[voisin2].x0)/6. / vol;
	double c_voisin2 = cross_product(Vector_3(F->centre,solide[part_1].x0),Vector_3(F->centre,solide[part_2].x0))*Vector_3(F->centre,solide[voisin1].x0)/6. / vol;
	if(vol > pow(10., -10.) && c_part_1 > 0. && c_part_2 > 0. && c_voisin1 > 0. && c_voisin2 > 0.) { //Stockage des particules du tetra et des coords bary si ok
	  (F->voisins).push_back(voisin1);
	  (F->voisins).push_back(voisin2);
	  (F->c_voisins).push_back(c_part_1);
	  (F->c_voisins).push_back(c_part_2);
	  (F->c_voisins).push_back(c_voisin1);
	  (F->c_voisins).push_back(c_voisin2);
	  tetra_ok = true;
	  break;
	}
      }
    }
    if(not(tetra_ok) && not(part_2 == -1)) { //Si marche pas avec part_1, on teste avec part_2
      int voisin1 = -1, voisin2 = -1; //Vont être les 2 autres particules composant le tetra associé à la face 
      //Boucle dans les faces de la première particule
      for(std::vector<int>::iterator G=solide[part_2].faces.begin();G!=solide[part_2].faces.end();G++) {
	if(not(faces[*G] == *F) && faces[*G].voisins[0] != part_2 && faces[*G].voisins[0] != -1 && voisin1 == -1)
	  voisin1 = faces[*G].voisins[0];
	else if(not(faces[*G] == *F) && faces[*G].voisins[1] != -1 && faces[*G].voisins[1] != part_2 && voisin1 == -1)
	  voisin1 = faces[*G].voisins[1];
      }
      for(std::vector<int>::iterator G=solide[part_2].faces.begin();G!=solide[part_2].faces.end();G++) {
	if(not(faces[*G] == *F) && faces[*G].voisins[0] != -1 && faces[*G].voisins[0] != part_1 && voisin1 != -1)
	  voisin2 = faces[*G].voisins[0];
	else if(not(faces[*G] == *F) && faces[*G].voisins[1] != -1 && faces[*G].voisins[1] != part_1 && voisin1 != -1)
	  voisin2 = faces[*G].voisins[1]; 

	//Vérifie que le centre de la face est dans le tetra et tetra pas aplati
	double vol = cross_product(Vector_3(solide[part_2].x0,solide[part_1].x0),Vector_3(solide[part_2].x0,solide[voisin1].x0))*Vector_3(solide[part_2].x0,solide[voisin2].x0)/6.; //Volume du tetra associé à la face
	//cout << F->id << endl;
	double c_part_1 = cross_product(Vector_3(F->centre,solide[part_2].x0),Vector_3(F->centre,solide[voisin1].x0))*Vector_3(F->centre,solide[voisin2].x0)/6. / vol;
	double c_part_2 = cross_product(Vector_3(F->centre,solide[part_1].x0),Vector_3(F->centre,solide[voisin1].x0))*Vector_3(F->centre,solide[voisin2].x0)/6. / vol;
	double c_voisin1 = cross_product(Vector_3(F->centre,solide[part_1].x0),Vector_3(F->centre,solide[part_2].x0))*Vector_3(F->centre,solide[voisin2].x0)/6. / vol;
	double c_voisin2 = cross_product(Vector_3(F->centre,solide[part_1].x0),Vector_3(F->centre,solide[part_2].x0))*Vector_3(F->centre,solide[voisin1].x0)/6. / vol;
        if(vol > pow(10., -10.) && c_part_1 > 0. && c_part_2 > 0. && c_voisin1 > 0. && c_voisin2 > 0.) { //Stockage des particules du tetra et des coords bary si ok
	  (F->voisins).push_back(voisin1);
	  (F->voisins).push_back(voisin2);
	  (F->c_voisins).push_back(c_part_1);
	  (F->c_voisins).push_back(c_part_2);
	  (F->c_voisins).push_back(c_voisin1);
	  (F->c_voisins).push_back(c_voisin2);
	  tetra_ok = true;
	  break;
	}
      }
    }
  }
}


/*std::vector<int> Solide::vertex_face(const int& particule, const int& voisin) {
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
}*/

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
  Forces_internes(dt);
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++)
    (P->second).Fi = (P->second).Fi + Forces_externes(t,T); //Reprendre calcul des forces externes
}

void Solide::stresses(){ //Calcul de la contrainte dans toutes les particules
  for(int i=0; i<faces.size(); i++){ //Calcul de la reconstruction sur chaque face
    faces[i].I_Dx = Vector_3(0., 0., 0.); //Remise à zéro. Si particule sur le bord, on a bien I_Dx = (0., 0., 0.)
    if(faces[i].BC != 1) {
       //not(faces[i].voisins[0] == -1) && not(faces[i].voisins[1] == -1)) { //Cad particule dans le bulk et pas sur le bord
      cout << i << endl;
      for(int j=0; i<faces[i].voisins.size() ; j++) {
	cout << "Pb ici ???" <<faces[i].c_voisins.size() << endl;
	cout << solide[faces[i].voisins[j]].Dx<< endl;
	faces[i].I_Dx = faces[i].I_Dx + faces[i].c_voisins[j] * solide[faces[i].voisins[j]].Dx;
      }
    }
  }
  
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    (P->second).discrete_gradient.col2 = Vector_3(0., 0., 0.);
    (P->second).discrete_gradient.col3 = Vector_3(0., 0., 0.);
    for(int i=0 ; i < (P->second).faces.size() ; i++){
      int f = (P->second).faces[i];
      Vector_3 nIJ = faces[i].normale;
      Matrix Dij_n(tens_sym(faces[f].I_Dx - (P->second).Dx,  nIJ) );
      (P->second).discrete_gradient += faces[i].S /  (P->second).V * Dij_n;
    }
    (P->second).contrainte = lambda * ((P->second).discrete_gradient - (P->second).epsilon_p).tr() * unit() + 2*mu * ((P->second).discrete_gradient - (P->second).epsilon_p);
	
    (P->second).seuil_elas = A; // + B * pow((P->second).def_plas_cumulee, n);

    if(((P->second).contrainte - H * (P->second).epsilon_p).VM() > (P->second).seuil_elas) { //On sort du domaine élastique.
      Matrix n_elas( 1. / (((P->second).contrainte).dev()).norme() * ((P->second).contrainte).dev() ); //Normale au domaine élastique de Von Mises
      double delta_p = (((P->second).contrainte - H * (P->second).epsilon_p).VM() - A) / (2*mu + H);
      //(P->second).def_plas_cumulee += delta_p;
      //(P->second).epsilon_p += delta_p * n_elas;
    }
  }

  /*if(plastifie)
    cout << "Plastification dans ce pas de temps !" << endl;*/
}


void Solide::Forces_internes(const double& dt){ //Calcul des forces pour chaque particule
  cout << "Avant calcul efforts !" << endl;
  stresses();
  for(std::map<int, Particule>::iterator P=solide.begin(); P!=solide.end(); P++){
    (P->second).Fi = Vector_3(0.,0.,0.);
    for(int i=0 ; i<(P->second).faces.size() ; i++){
      if(faces[i].BC != 0){ //On prend pas les faces au bord car il n'y a pas de forces internes dedans
	int num_face = (P->second).faces[i]; //Numéro de la face dans l'ensemble des faces contenu dans le solide
	int part_1 = faces[num_face].voisins[0];
	double c_part_1 = faces[num_face].c_voisins[0];
	int part_2 = faces[num_face].voisins[1];
	double c_part_2 = faces[num_face].c_voisins[1];
	int aux_1 = faces[num_face].voisins[2];
	double c_aux_1 = faces[num_face].c_voisins[2];
	int aux_2 = faces[num_face].voisins[3];
	double c_aux_2 = faces[num_face].c_voisins[3];
	//Sortir le sens de toutes les forces comme il faut...
	Vector_3 nIJ = faces[num_face].normale;
	//Ajouter les directions des forces avec particules aux_1 et aux_2
	Vector_3 n_aux_1 = Vector_3((P->second).x0, solide[aux_1].x0) / sqrt(Vector_3((P->second).x0, solide[aux_1].x0).squared_length());
	Vector_3 n_aux_2 = Vector_3((P->second).x0, solide[aux_2].x0) / sqrt(Vector_3((P->second).x0, solide[aux_2].x0).squared_length());
	(P->second).Fi = (P->second).Fi + faces[num_face].S * (c_part_1 * solide[part_1].contrainte + c_part_2 * solide[part_2].contrainte) * nIJ + faces[num_face].S * c_aux_1 * solide[aux_1].contrainte * n_aux_1 + faces[num_face].S * c_aux_2 * solide[aux_2].contrainte * n_aux_2;
      }
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

  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    Ep += 0.5 * contraction_double((P->second).contrainte, (P->second).discrete_gradient - (P->second).epsilon_p) * (P->second).V; //+ B * pow(((P->second)).def_plas_cumulee, 1. + n) / (n + 1.) + A * ((P->second)).def_plas_cumulee ) * ((P->second)).V;
  }
  //cout << "Energie potentielle : " << Ep << endl;
  return Ep;
}

//Reprendre le calcul de la CFL !!!
double Solide::pas_temps(const double& t, const double& T, const double& cfls, const double& E, const double& nu, const double& rhos){
  double eps = 1e-14;//std::numeric_limits<double>::epsilon();
  double dt = std::numeric_limits<double>::infinity();
  //Restriction CFL liee aux forces internes
  double cs = sqrt(E*(1.-nu)/rhos/(1.+nu)/(1.-2.*nu));
  //Calcul du rayon de la sphère inscrite
  double sigma = 100000.;
  for(int i=0;i<size();i++){
    for(int j=0;j<solide[i].faces.size();j++){
      sigma = min(sigma,faces[solide[i].faces[j]].D0);
    }
  }
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(int j=0;j<(P->second).faces.size();j++){
      if(faces[(P->second).faces[j]].voisins[0] >=0 && faces[(P->second).faces[j]].voisins[1] >= 0){
	dt = min(dt,cfls*faces[(P->second).faces[j]].D0/cs);
      }
    }
  }
  dt = min(dt,T-t);
  return dt;
}

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
    for(std::vector<int>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      vtk << faces[*F].vertex.size();
      for(int k=0 ; k<faces[*F].vertex.size() ; k++)
	vtk << " " << compteur_vertex + (faces[*F].vertex)[k];
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
    for(std::vector<int>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      vtk << (P->second).Dx << endl;
  }
  vtk << "\n";
  //Vitesse
  vtk << "VECTORS vitesse double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      vtk << (P->second).u << endl;
  }
  vtk << "\n";
  //Contrainte
  vtk << "TENSORS contraintes double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
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
    for(std::vector<int>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
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
    for(std::vector<int>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
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
    for(std::vector<int>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      vtk << (P->second).def_plas_cumulee << endl;
  }
  vtk << "\n";
  vtk.close();
}

#endif
