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
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    solide[P->id]= *P;
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
    int id; //Numéro du vertex
    double x,y,z;
    stm >> id >> x >> y >> z;
    vertex.push_back(Vertex(Point_3(x,y,z), id)); //Vertex sont donnés dans l'ordre
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
    p.vertices.push_back(vertex[v1].pos);
    p.vertices.push_back(vertex[v2].pos);
    p.vertices.push_back(vertex[v3].pos);
    p.vertices.push_back(vertex[v4].pos);

    //Calcul des quantités volumiques (liées particule)
    p.barycentre(); //Calcul du barycentre
    p.volume(); //calcul du volume
    p.m = rho * p.V;

    //Ajout de la particule dans le solide
    solide.push_back(p);
  }

  //cout << solide.begin()->first <<" " << solide.end()->first << endl;
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
    int BC; //Condition de bord
    int part_1, part_2; //Le numéro des particules voisines
    stm >> id >> v1 >> v2 >> v3 >> BC >> part_1 >> part_2; //Numéro de la face, des vertex + 1 et des voisins (bon numéro)
    Face f;
    f.id = id;
    f.vertex.push_back(v1); //Ajout du numéro des vertex
    f.vertex.push_back(v2);
    f.vertex.push_back(v3);
    f.BC = BC;
    /*if(part_1 == -1 || part_2 == -1)
      f.BC = -1; //face au bord
    else {
      f.BC = 0; //Face pas au bord
      f.D0 = sqrt(Vector_3(solide[part_1].x0, solide[part_2].x0).squared_length());
      }*/
    if(part_1 >=0) { //Ajout du numéro des voisins dans la face
      f.voisins.push_back(part_1);
      f.voisins.push_back(part_2);
    }
    else { //La première valeur de voisin est la seule particule qui contient la face sur le bord
      f.voisins.push_back(part_2);
      f.voisins.push_back(part_1);
    }
    //bool calcul_normales = false;
    if(part_1 >= 0) {
      solide[part_1].faces.push_back(f.id); //Ajout du numéro de la face dans la liste ds voisins de chaque particule
      if(f.BC != 0) {
	if(f.BC == 1)
	  solide[part_1].BC = f.BC;
	else if(solide[part_1].BC <= 0)
	  solide[part_1].BC = f.BC;
      }
    }
    if(part_2 >= 0) {
      solide[part_2].faces.push_back(f.id); //Ajout du numéro de la face dans la liste ds voisins de chaque particule
      if(f.BC != 0) {
	if(f.BC == 1)
	  solide[part_1].BC = f.BC;
	else if(solide[part_1].BC <= 0)
	  solide[part_1].BC = f.BC;
      }
    }    
    f.comp_quantities(vertex[v1].pos, vertex[v2].pos, vertex[v3].pos); //Calcul de la normale sortante, surface et barycentre face

    //Vérification du sens de la normale
    if(part_1 != -1 && part_2 != -1) { //Face pas au bord
      if(Vector_3(solide[part_1].x0, solide[part_2].x0) * f.normale  < 0.)
	f.normale = -f.normale;
    }
    if(part_1 == -1 || part_2 == -1) { //Face au bord
      Vector_3 bonne_direction = Vector_3(solide[f.voisins[0]].x0, f.centre);
      if(bonne_direction * f.normale < 0.)
	f.normale = -f.normale;
    }
    faces.push_back(f);
  }

  //Calcul du tetrahèdre associé à chaque face pour le calcul du gradient
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){ //Boucle sur toutes les faces
    if(F->BC == 0) {
      //cout << "Face : " << F->id << endl;
      bool test = voisins_face(F->id);
      if(not(test)) {
	//cout << "Face : " << F->id << " Pas de tetra associe a une face" << endl;
	throw std::invalid_argument( "Pas de tetra associe a une face" );
	/*cout << "Centre Face : " << F->centre << endl;
	cout << "Barycentre Voisin A : " << solide[F->voisins[0]].x0 << endl;
	cout << "Barycentre Voisin A : " << solide[F->voisins[0]].x0 << endl;
	cout << "Barycentre Voisin B : " << solide[F->voisins[1]].x0 << endl;*/
      }
    }
  }
}

void Solide::Init(const char* s, const bool& rep, const int& numrep, const double& rho){ //Pour gmsh

}

bool Solide::voisins_face(int num_face) {
  int part_1 = faces[num_face].voisins[0];
  int part_2 = faces[num_face].voisins[1];
  //cout << part_1 << " " << part_2 << endl;

  vector<int> tous_voisins; //Va être rempli des voisins des 2 particules qui peuvent être candidats pour former le tetra associé à la face !

  for(std::vector<int>::iterator G=solide[part_1].faces.begin();G!=solide[part_1].faces.end();G++){
    //cout << "Num face teste : " << *G << endl;
    if(not(faces[*G] == faces[num_face]) && faces[*G].voisins[0] != part_1 && faces[*G].voisins[0] != -1)
      tous_voisins.push_back(faces[*G].voisins[0]);
    else if(not(faces[*G] == faces[num_face]) && faces[*G].voisins[1] != part_1 && faces[*G].voisins[1] != -1)
      tous_voisins.push_back(faces[*G].voisins[1]);
  }
  //Boucle sur deuxième groupe de particules
  for(std::vector<int>::iterator G=solide[part_2].faces.begin();G!=solide[part_2].faces.end();G++){
    //cout << "Num face teste : " << *G << endl;
    if(not(faces[*G] == faces[num_face]) && faces[*G].voisins[0] != part_2 && faces[*G].voisins[0] != -1)
      tous_voisins.push_back(faces[*G].voisins[0]);
    else if(not(faces[*G] == faces[num_face]) && faces[*G].voisins[1] != part_2 && faces[*G].voisins[1] != -1)
      tous_voisins.push_back(faces[*G].voisins[1]);
  }
  //Tous_voisins rempli a priori

  bool tetra_ok = false;
  for(std::vector<int>::iterator G=tous_voisins.begin();G!=tous_voisins.end()-1;G++){
    if(tetra_ok)
      break;
    for(std::vector<int>::iterator I=G + 1;I!=tous_voisins.end();I++){
      int voisin1 = *G;
      int voisin2 = *I;
      double vol = abs(cross_product(Vector_3(solide[part_1].x0,solide[part_2].x0),Vector_3(solide[part_1].x0,solide[voisin1].x0))*Vector_3(solide[part_1].x0,solide[voisin2].x0)/6.); //Volume du tetra associé à la face
      if(vol > pow(10., -8.)) {
	double c_part_1 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_2].x0, solide[part_1].x0) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0) ));
	double c_part_2 = (Vector_3(solide[part_1].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_1].x0, solide[part_2].x0) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0) ));
	double c_voisin1 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_2].x0, solide[part_1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_2].x0, solide[voisin1].x0) * cross_product(Vector_3(solide[part_2].x0, solide[part_1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0) ));
	double c_voisin2 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[part_1].x0)) ) / (Vector_3(solide[part_2].x0, solide[voisin2].x0) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[part_1].x0) ));

	faces[num_face].voisins.push_back(voisin1);
	faces[num_face].voisins.push_back(voisin2);
	faces[num_face].c_voisins.push_back(c_part_1);
	faces[num_face].c_voisins.push_back(c_part_2);
	faces[num_face].c_voisins.push_back(c_voisin1);
	faces[num_face].c_voisins.push_back(c_voisin2);
	tetra_ok = true;
	break;	
      }
    }
  }
  return tetra_ok;
}

void Solide::Solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    P->solve_position(dt, flag_2d, t, T);
  }
}

void Solide::Solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    P->solve_vitesse(dt, flag_2d, Amort, t , T);
  }
}

void Solide::Forces(const int& N_dim, const double& dt, const double& t, const double& T){
  Forces_internes(dt, t);
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++)
    P->Fi = P->Fi; // + Forces_externes(t,T); //Reprendre calcul des forces externes
}

void Solide::stresses(const double& t){ //Calcul de la contrainte dans toutes les particules
  for(int i=0; i<faces.size(); i++){ //Calcul de la reconstruction sur chaque face
    if(t > 0.)
      faces[i].I_Dx = Vector_3(0., 0., 0.); //Remise à zéro. Si particule sur le bord, on a bien I_Dx = (0., 0., 0.)
    //cout << "BC : " << faces[i].BC << endl;
    //Vector_3 test_pos(0., 0., 0.);
    if(faces[i].BC == 1) {
      if(t > 0.)
	faces[i].I_Dx = solide[faces[i].voisins[0]].Dx; //Dirichlet BC imposée fortement dans Mka ! old...
      //cout << faces[i].I_Dx.vec[2] << endl;
      //faces[i].I_Dx = displacement_BC(faces[i].centre, solide[faces[i].voisins[0]].Dx, t, 0.);
      //if(t < pow(10., -8.))
      //faces[i].I_Dx.vec[2] = displacement_BC_bis(faces[i].centre, solide[faces[i].voisins[0]].Dx, t, 0.); //BC de Dirichlet
    }
    else if(faces[i].BC == -1) {
      faces[i].I_Dx = solide[faces[i].voisins[0]].Dx; //Dirichlet BC imposée fortement dans Mka ! old...
      //faces[i].I_Dx = displacement_BC(faces[i].centre, solide[faces[i].voisins[0]].Dx, t, 0.);
      //faces[i].I_Dx.vec[2] = displacement_BC_bis(faces[i].centre, solide[faces[i].voisins[0]].Dx, t, 0.);
    }
    else if(faces[i].BC == 0) { //Cad particule dans le bulk
      for(int j=0; j<faces[i].voisins.size() ; j++) {
	faces[i].I_Dx = faces[i].I_Dx + faces[i].c_voisins[j] * solide[faces[i].voisins[j]].Dx;
	//test_pos = test_pos + faces[i].c_voisins[j] * Vector_3(Point_3(0., 0., 0.), solide[faces[i].voisins[j]].x0);
      }
      //faces[i].I_Dx = faces[i].c_voisins[0] * solide[faces[i].voisins[0]].Dx + faces[i].c_voisins[1] * solide[faces[i].voisins[1]].Dx;
      //faces[i].I_Dx = (solide[faces[i].voisins[0]].Dx + solide[faces[i].voisins[1]].Dx) / 2.;

      /*if(abs(faces[i].I_Dx[2] - 4. / 3. * faces[i].centre.z()) > 0.00001)
	cout << "Problème reconstruction sur face : " << i << endl;*/
      /*if(sqrt((test_pos - Vector_3(Point_3(0.,0.,0.),faces[i].centre)).squared_length()) > pow(10., -11.))
	cout << "Problème reconstruction barycentre face : " << i << endl;*/
    }
  }
  
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    P->discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    P->discrete_gradient.col2 = Vector_3(0., 0., 0.);
    P->discrete_gradient.col3 = Vector_3(0., 0., 0.);
    //Matrix test;
    //Vector_3 test_vec;
    //Test
    /*if(abs(P->Dx[2] - 4. / 3. * P->x0.z()) > pow(10., -5.))
      cout << "Problème reconstruction sur cellule : " << P->id << endl;*/
    for(int i=0 ; i < P->faces.size() ; i++){
      int f = P->faces[i];
      Vector_3 nIJ = faces[f].normale;
      if(faces[f].BC == 0 && nIJ * Vector_3(P->x0, faces[f].centre) < 0.)
	  nIJ = -nIJ; //Normale pas dans le bon sens...
      /*Vector_3 nIJ = Vector_3(P->x0, faces[f].centre);
	nIJ = nIJ / sqrt(nIJ.squared_length());*/
      /*if(abs(nIJ.squared_length() - 1.) > pow(10., -10.))
	cout << "Pas la bonne norme !!!!" << endl;*/
      if(faces[f].BC >= 0) { //Car conditions de Neumann homogènes sur bord en -1
	/*int voisin;
	if(P->id == faces[f].voisins[0])
	  voisin = faces[f].voisins[1];
	else if(P->id == faces[f].voisins[1])
	voisin = faces[f].voisins[0];*/
	//Matrix Dij_n( tens_sym(solide[voisin].Dx - P->Dx,  nIJ) / 2. ); //OK Voronoi
	Matrix Dij_n(tens_sym(faces[f].I_Dx - P->Dx,  nIJ) ); //Tetra
	P->discrete_gradient += faces[f].S /  P->V * Dij_n;
      }
      //P->discrete_gradient += tens_sym(solide[voisin].Dx - P->Dx, nIJ);
      //test = test + faces[f].S /  P->V * tens(Vector_3(P->x0,faces[f].centre),  nIJ);
      //test_vec = test_vec + faces[f].S * nIJ;
      /*if(P->faces.size() != 4)
      cout << "Nbr faces : " << P->faces.size() << endl;
      if(faces[f].S<0.){
      cout << "S=" << faces[f].S << endl;*/
    }
      /*if(P->discrete_gradient.col3[2] > pow(10., -6.))
	cout << "Num element : " << P->id << " def : " << P->discrete_gradient.col3[2] << endl;*/

    /*if(test_vec.squared_length() > pow(10., -5.)) {
    //cout << "MTF : " << test_vec << endl;
    cout << "Num element : " << P->id << endl;
      for(int i=0 ; i < P->faces.size() ; i++){
	cout << vertex[faces[i].vertex[0]].pos << endl;
	cout << vertex[faces[i].vertex[1]].pos << endl;
	cout << vertex[faces[i].vertex[2]].pos << endl;
	cout << endl;
	}*/
    
/*if(sqrt(contraction_double(test - Matrix(Vector_3(1.,0.,0.), Vector_3(0.,1.,0.), Vector_3(0.,0.,1.)), test - Matrix(Vector_3(1.,0.,0.), Vector_3(0.,1.,0.), Vector_3(0.,0.,1.)))) > pow(10.,-5.))
      cout << "Problème sur tenseur identité !" << endl;
    cout << test.col1 << endl;
    cout << test.col2 << endl;
    cout << test.col3 << endl;
    cout << "V=" << P->V <<  endl;*/
    P->contrainte = lambda * (P->discrete_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_gradient - P->epsilon_p);
    P->seuil_elas = A; // + B * pow(P->def_plas_cumulee, n);

    if((P->contrainte - H * P->epsilon_p).VM() > P->seuil_elas) { //On sort du domaine élastique.
      Matrix n_elas( 1. / ((P->contrainte).dev()).norme() * (P->contrainte).dev() ); //Normale au domaine élastique de Von Mises
      double delta_p = ((P->contrainte - H * P->epsilon_p).VM() - A) / (2*mu + H);
      //P->def_plas_cumulee += delta_p;
      //P->epsilon_p += delta_p * n_elas;
    }
  }
}


void Solide::Forces_internes(const double& dt, const double& t){ //Calcul des forces pour chaque particule
  stresses(t);
  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++) //Remet à zéro toutes les forces
    P->Fi = Vector_3(0.,0.,0.);
  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++){
    for(int i=0 ; i<P->faces.size() ; i++){
      int num_face = P->faces[i]; //Numéro de la face dans l'ensemble des faces contenu dans le solide
      int part_1 = faces[num_face].voisins[0];
      int part_2 = faces[num_face].voisins[1];
      if(faces[num_face].BC == 0) { // && not(part_1 == -1 || part_2 == -1)){ //On prend pas les faces au bord car il n'y a pas de forces internes dedans
	double c_part_1 = faces[num_face].c_voisins[0];
	double c_part_2 = faces[num_face].c_voisins[1];
	int aux_1 = faces[num_face].voisins[2];
	double c_aux_1 = faces[num_face].c_voisins[2];
	int aux_2 = faces[num_face].voisins[3];
	double c_aux_2 = faces[num_face].c_voisins[3];
	//cout << "coords bary : " <<c_part_1 << " " << c_part_2 << " " << c_aux_1 << " " << c_aux_2 << endl;
	
	//Sortir le sens de toutes les forces comme il faut...
	Vector_3 nIJ = faces[num_face].normale;
	int voisin;
	if(P->id == faces[num_face].voisins[0])
	  voisin = faces[num_face].voisins[1];
	else if(P->id == faces[num_face].voisins[1])
	  voisin = faces[num_face].voisins[0];
	if(nIJ * Vector_3(P->x0, faces[num_face].centre) < 0.)
	  nIJ = -nIJ; //Normale pas dans le bon sens...

	//P->Fi = P->Fi +  2 * mu * (solide[voisin].Dx - P->Dx); //OK flux à 2 points
	//P->Fi = P->Fi + faces[num_face].S * (solide[part_1].contrainte + solide[part_2].contrainte) / 2. * nIJ; //OK Voronoi
	P->Fi = P->Fi + faces[num_face].S * (c_part_2 * solide[part_1].contrainte + c_part_1 * solide[part_2].contrainte) * nIJ + faces[num_face].S * c_aux_1 * P->contrainte * nIJ + faces[num_face].S * c_aux_2 *P->contrainte * nIJ;
	solide[aux_1].Fi = solide[aux_1].Fi - faces[num_face].S * c_aux_1 * P->contrainte * nIJ;
	solide[aux_2].Fi = solide[aux_2].Fi - faces[num_face].S * c_aux_2 * P->contrainte * nIJ;
      }
      /*else if(faces[num_face].BC == 1 && t < pow(10., -8.)) { //Pas de DDL sur face avec BC de Neuman homogène et Dirichlet que à t=0. Ensuite Neuman homogène
	int part = faces[num_face].voisins[0];
	Vector_3 nIJ = faces[num_face].normale;
	P->Fi = P->Fi + faces[num_face].S * pow(10., 7.) * nIJ;//solide[part].contrainte * nIJ;
	}*/
    }
    /*cout << "Particule :" << P->first << endl;
    cout << "Force : " << P->Fi << endl; */
  }
}

const double Solide::Energie(){
  return Energie_cinetique()+Energie_potentielle();
}

const double Solide::Energie_cinetique(){
  double E = 0.;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    E += 1./2. * P->m * (P->u + P->u_prev) * (P->u + P->u_prev) / 4.;
    //E += 1./2. * P->m * P->u * P->u;
  }
  //cout << "Energie cinetique : " << E << endl;
  return E;
  //return 0.;
}

const double Solide::Energie_potentielle(){
  double Ep = 0.;

  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    //Juste flux à 2 points
    /*for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++){
      int voisin;
      int num_face = *F;
      if(faces[num_face].BC == 0) {
	if(P->id == faces[num_face].voisins[0])
	  voisin = faces[num_face].voisins[1];
	else if(P->id == faces[num_face].voisins[1])
	  voisin = faces[num_face].voisins[0];
	Ep += 0.25 * 2. * mu * (solide[voisin].Dx - P->Dx) * (solide[voisin].Dx - P->Dx);
      }
      }*/
    //Ep += 0.5 * contraction_double(P->contrainte, P->discrete_gradient - P->epsilon_p) * P->V; //OK Voronoi
    Ep += 0.5 * contraction_double(P->contrainte, P->discrete_gradient - P->epsilon_p) * P->V; /*+ B * pow(((P->second)).def_plas_cumulee, 1. + n) / (n + 1.)*/ + A * P->def_plas_cumulee * P->V;
    //}
  }
  //cout << "Energie potentielle : " << Ep << endl;
  //return 0.;
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
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(int j=0;j<P->faces.size();j++){
      if(faces[P->faces[j]].voisins[0] >=0 && faces[P->faces[j]].voisins[1] >= 0){
	dt = min(dt,cfls*faces[P->faces[j]].D0/cs);
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
  int nb_points = 3 * 4 * nb_part;
  int nb_faces = 4 * nb_part;
  int size = 3 * nb_faces;
  /*for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      size += F->nb_vertex; //Egale à nb_points au final ?
  }*/
  size += nb_faces; //Pk ? Mystère...
    
  //Initialisation du fichier vtk
  vtk << "# vtk DataFile Version 3.0" << endl;
  vtk << "#Simulation Euler" << endl;
  vtk << "ASCII" << endl;
  vtk<< "\n";
  vtk << "DATASET UNSTRUCTURED_GRID" << endl;
  vtk << "POINTS " << nb_points << " DOUBLE" << endl;
    
  //Sortie des points
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
      for(std::vector<int>::iterator V=faces[*F].vertex.begin();V!=faces[*F].vertex.end();V++)
	vtk << P->mvt_t(vertex[*V].pos) << endl;
    }
  }
  vtk << endl;
    
  //Sortie des faces
  int point_tmp=0;
  vtk << "CELLS " << nb_faces << " " << size << endl;
  int compteur_vertex = 0;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
      vtk << faces[*F].vertex.size(); //Vaut 3
      for(int k=0 ; k<faces[*F].vertex.size() ; k++)
	vtk << " " << compteur_vertex + k; //(faces[*F].vertex)[k];
      vtk << endl;
      compteur_vertex += faces[*F].vertex.size();
    }
    //vtk << endl;
  }
  vtk << endl;
  vtk << "CELL_TYPES " << nb_faces << endl;
  for(int i=0;i<nb_faces;i++){
    //vtk << 7 << endl;
    vtk << 5 << endl; //Pour triangle
  }
  vtk << "\n";
  vtk << "CELL_DATA " << nb_faces << endl;
  //Deplacement
  vtk << "VECTORS deplacement double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      vtk << P->Dx << endl;
  }
  vtk << "\n";
  //Vitesse
  vtk << "VECTORS vitesse double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      vtk << P->u << endl;
  }
  vtk << "\n";
  //Normale
  vtk << "VECTORS normale double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
    vtk << faces[*F].normale << endl;
  }
  vtk << "\n";
  //Contrainte
  vtk << "TENSORS contraintes double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
      vtk << P->contrainte.col1[0] << " " << P->contrainte.col1[1] << " " << P->contrainte.col1[2] << endl;
      vtk << P->contrainte.col2[0] << " " << P->contrainte.col2[1] << " " << P->contrainte.col2[2] << endl;
      vtk << P->contrainte.col3[0] << " " << P->contrainte.col3[1] << " " << P->contrainte.col3[2] << endl;
    }
  }
  vtk << "\n";
  //Déformations
  vtk << "TENSORS deformations double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
      vtk << P->discrete_gradient.col1[0] << " " << P->discrete_gradient.col1[1] << " " << P->discrete_gradient.col1[2] << endl;
      vtk << P->discrete_gradient.col2[0] << " " << P->discrete_gradient.col2[1] << " " << P->discrete_gradient.col2[2] << endl;
      vtk << P->discrete_gradient.col3[0] << " " << P->discrete_gradient.col3[1] << " " << P->discrete_gradient.col3[2] << endl;
    }
  }
  vtk << "\n";
  //Epsilon_p
  vtk << "TENSORS epsilon_p double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
      vtk << P->epsilon_p.col1[0] << " " << P->epsilon_p.col1[1] << " " << P->epsilon_p.col1[2] << endl;
      vtk << P->epsilon_p.col2[0] << " " << P->epsilon_p.col2[1] << " " << P->epsilon_p.col2[2] << endl;
      vtk << P->epsilon_p.col3[0] << " " << P->epsilon_p.col3[1] << " " << P->epsilon_p.col3[2] << endl;
    }
  }
  vtk << "\n";
  //Deformation plastique cumulée
  vtk << "SCALARS p double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      vtk << P->def_plas_cumulee << endl;
  }
  vtk << "\n";
  //Numéro des particules
  /*vtk << "SCALARS id double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      vtk << P->id << endl;
  }
  vtk << endl;*/
  vtk.close();
}

#endif
