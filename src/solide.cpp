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
//#include <eigen3/Eigen/LU> //Sert pour inversion systeème linéaire pour condition Neumann
#include <eigen3/Eigen/Dense>
#include <cmath>
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

void Solide::Init(const char* s1, const char* s2, const char* s3, const bool& rep, const int& numrep, const double& rho){ //pour Tetgen
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
    p.vertices.push_back(v1);
    p.vertices.push_back(v2);
    p.vertices.push_back(v3);
    p.vertices.push_back(v4);

    //Calcul des quantités volumiques (liées particule)
    p.barycentre(this, 4); //Calcul du barycentre
    p.volume(this, 4); //calcul du volume
    p.m = rho * p.V;

    //Ajout de la particule dans le solide
    solide.push_back(p);
  }
  
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
    f.type = 2; //Triangle pour tegen
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
    f.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face

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

void Solide::Init(const char* s1, const bool& rep, const int& numrep, const double& rho){ //Pour gmsh
  std::ifstream maillage(s1,ios::in); //Ouverture du maillage
  if(not(maillage))
    throw std::invalid_argument( "Ouverture du maillage ratee !" );

  //Importation des vertex
  string s;
  char aux;
  string ligne;
  while(getline(maillage, ligne) && ligne != "$Nodes") {} //On fait rien...
  getline(maillage, ligne); //Nombre de vertex. Pas utile.

  //Importation des Vertex
  while(getline(maillage, ligne) && ligne != "$EndNodes") {
    istringstream  stm(ligne);
    int id; //Numéro du vertex
    double x,y,z;
    stm >> id >> x >> y >> z;
    vertex.push_back(Vertex(Point_3(x,y,z), id-1)); //Vertex sont donnés dans l'ordre
  }

  //Importation des Particules
  int type; //Celui des particules
  getline(maillage, ligne); //$Elements
  getline(maillage, ligne); //Nbr Elements
  while(getline(maillage, ligne) && ligne != "$EndElements") {
    istringstream  stm(ligne);
    int id,nbr_tag,tag_1,tag_2;
    stm >> id >> type;

    //Importation des faces au bord pour avoir les BC !
    if(type == 2) { //Triangle donc sur bord
      int v1,v2,v3;
      stm >> nbr_tag >> tag_1 >> tag_2 >> v1 >> v2 >> v3; //tag_2 est la BC
      Face F;
      F.vertex.push_back(v1 - 1);
      F.vertex.push_back(v2 - 1);
      F.vertex.push_back(v3 - 1);
      F.type = 2;
      if(tag_2 == 11 || tag_2 == 33) //Dirichlet
	F.BC = 1;
      else if(tag_2 == 20 || tag_2 == 24 || tag_2 == 28 || tag_2 == 32) //Neumann
	F.BC = -1;
      /*if(tag_2 == 6 || tag_2 == 28) //Dirichlet
	F.BC = -1; //Tout en Neumann pour ce test
      else if(tag_2 == 15 || tag_2 == 19 || tag_2 == 23 || tag_2 == 27) //Neumann
      F.BC = -1;*/
      F.id = faces.size();
      //F.voisins.push_back(F.id); F.voisins.push_back(-1); //Car face au bord
      F.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
      if(F.normale * Vector_3(F.centre, vertex[0].pos) < 0.)
	  F.normale = -F.normale;
      faces.push_back(F);
    }
    else if(type == 4) { //Tetra
      int v1,v2,v3,v4;
      stm >> nbr_tag >> tag_1 >> tag_2 >> v1 >> v2 >> v3 >> v4;
      //Ajout des vertex de la particule
      Particule p;
      p.id = solide.size();
      p.vertices.push_back(v1 - 1);
      p.vertices.push_back(v2 - 1);
      p.vertices.push_back(v3 - 1);
      p.vertices.push_back(v4 - 1);

      //Calcul des quantités volumiques (liées particule)
      p.barycentre(this, type); //Calcul du barycentre
      p.volume(this, type); //calcul du volume
      p.m = rho * p.V;

      //Ajout de la particule dans le solide
      solide.push_back(p);
    }
    else if(type == 5) { //Hexa
      int v1,v2,v3,v4,v5,v6,v7,v8;
      stm >> nbr_tag >> tag_1 >> tag_2 >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7 >> v8;
      //Ajout des vertex de la particule
      Particule p;
      p.id = solide.size();
      p.vertices.push_back(v1 - 1);
      p.vertices.push_back(v2 - 1);
      p.vertices.push_back(v3 - 1);
      p.vertices.push_back(v4 - 1);
      p.vertices.push_back(v5 - 1);
      p.vertices.push_back(v6 - 1);
      p.vertices.push_back(v7 - 1);
      p.vertices.push_back(v8 - 1);

      //Calcul des quantités volumiques (liées particule)
      p.barycentre(this, type); //Calcul du barycentre
      p.volume(this, type); //calcul du volume
      p.m = rho * p.V;

      //Ajout de la particule dans le solide
      solide.push_back(p);
    }
  }

  //Création des faces
  if(type == 5) { //Hexa
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      //face 1
      Face face1;
      face1.vertex.push_back(P->vertices[0]);
      face1.vertex.push_back(P->vertices[3]);
      face1.vertex.push_back(P->vertices[2]);
      face1.vertex.push_back(P->vertices[1]);
      face1.type = 3;
      if(not(face_existe(face1))) { //Ajout de la face dans l'ensemble des faces du Solide
	face1.id = faces.size();
	face1.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face1.normale * Vector_3(face1.centre, vertex[P->vertices[1]].pos) < 0.)
	  face1.normale = -face1.normale;
	faces.push_back(face1);
      }

      //face 2
      Face face2;
      face2.vertex.push_back(P->vertices[0]);
      face2.vertex.push_back(P->vertices[4]);
      face2.vertex.push_back(P->vertices[7]);
      face2.vertex.push_back(P->vertices[3]);
      face2.type = 3;
      if(not(face_existe(face2))) { //Ajout de la face dans l'ensemble des faces du Solide
	face2.id = faces.size();
	face2.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face2.normale * Vector_3(face2.centre, vertex[P->vertices[0]].pos) < 0.)
	  face2.normale = -face2.normale;
	faces.push_back(face2);
      }

      //face 3
      Face face3;
      face3.vertex.push_back(P->vertices[4]);
      face3.vertex.push_back(P->vertices[5]);
      face3.vertex.push_back(P->vertices[6]);
      face3.vertex.push_back(P->vertices[7]);
      face3.type = 3;
      if(not(face_existe(face3))) { //Ajout de la face dans l'ensemble des faces du Solide
	face3.id = faces.size();
	face3.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face3.normale * Vector_3(face3.centre, vertex[P->vertices[5]].pos) < 0.)
	  face3.normale = -face3.normale;
	faces.push_back(face3);
      }

      //face 4
      Face face4;
      face4.vertex.push_back(P->vertices[5]);
      face4.vertex.push_back(P->vertices[1]);
      face4.vertex.push_back(P->vertices[2]);
      face4.vertex.push_back(P->vertices[6]);
      face4.type = 3;
      if(not(face_existe(face4))) { //Ajout de la face dans l'ensemble des faces du Solide
	face4.id = faces.size();
	face4.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face4.normale * Vector_3(face4.centre, vertex[P->vertices[2]].pos) < 0.)
	  face4.normale = -face4.normale;
	faces.push_back(face4);
      }

      //face 5
      Face face5;
      face5.vertex.push_back(P->vertices[7]);
      face5.vertex.push_back(P->vertices[6]);
      face5.vertex.push_back(P->vertices[2]);
      face5.vertex.push_back(P->vertices[3]);
      face5.type = 3;
      if(not(face_existe(face5))) { //Ajout de la face dans l'ensemble des faces du Solide
	face5.id = faces.size();
	face5.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face5.normale * Vector_3(face5.centre, vertex[P->vertices[2]].pos) < 0.)
	  face5.normale = -face5.normale;
	faces.push_back(face5);
      }

      //face 6
      Face face6;
      face6.vertex.push_back(P->vertices[0]);
      face6.vertex.push_back(P->vertices[1]);
      face6.vertex.push_back(P->vertices[5]);
      face6.vertex.push_back(P->vertices[4]);
      face6.type = 3;
      if(not(face_existe(face6))) { //Ajout de la face dans l'ensemble des faces du Solide
	face6.id = faces.size();
	face6.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face6.normale * Vector_3(face6.centre, vertex[P->vertices[4]].pos) < 0.)
	  face6.normale = -face6.normale;
	faces.push_back(face6);
      }
    }
  }
  else if(type == 4) { //Pour Tetra
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      //face 1
      Face face1;
      face1.vertex.push_back(P->vertices[0]);
      face1.vertex.push_back(P->vertices[3]);
      face1.vertex.push_back(P->vertices[2]);
      face1.type = 2;
      face1.BC = 0;
      if(not(face_existe(face1))) { //Ajout de la face dans l'ensemble des faces du Solide
	face1.id = faces.size();
	face1.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face1.normale * Vector_3(face1.centre, vertex[P->vertices[0]].pos) < 0.)
	  face1.normale = -face1.normale;
	faces.push_back(face1);
      }

      //face 2
      Face face2;
      face2.vertex.push_back(P->vertices[3]);
      face2.vertex.push_back(P->vertices[1]);
      face2.vertex.push_back(P->vertices[2]);
      face2.BC = 0;
      face2.type = 2;
      if(not(face_existe(face2))) { //Ajout de la face dans l'ensemble des faces du Solide
	face2.id = faces.size();
	face2.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face2.normale * Vector_3(face2.centre, vertex[P->vertices[3]].pos) < 0.)
	  face2.normale = -face2.normale;
	faces.push_back(face2);
      }

      //face 3
      Face face3;
      face3.vertex.push_back(P->vertices[1]);
      face3.vertex.push_back(P->vertices[0]);
      face3.vertex.push_back(P->vertices[2]);
      face3.BC = 0;
      face3.type = 2;
      if(not(face_existe(face3))) { //Ajout de la face dans l'ensemble des faces du Solide
	face3.id = faces.size();
	face3.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face3.normale * Vector_3(face3.centre, vertex[P->vertices[1]].pos) < 0.)
	  face3.normale = -face3.normale;
	faces.push_back(face3);
      }

      //face 4
      Face face4;
      face4.vertex.push_back(P->vertices[0]);
      face4.vertex.push_back(P->vertices[1]);
      face4.vertex.push_back(P->vertices[3]);
      face4.BC = 0;
      face4.type = 2;
      if(not(face_existe(face4))) { //Ajout de la face dans l'ensemble des faces du Solide
	face4.id = faces.size();
	face4.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
	if(face4.normale * Vector_3(face4.centre, vertex[P->vertices[0]].pos) < 0.)
	  face4.normale = -face4.normale;
	faces.push_back(face4);
      }
    }
  }

  //cout << "nombre particules : " << solide.size() << endl;
  //cout << "Nombre total de faces : " << faces.size() << endl;

  //Création des connectivités entre éléments
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){ //Boucle sur toutes les faces
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      if(P->contient_face(*F)) { //On ajoute les numéros de la face dans la particule et réciproquement
	//cout << "Particule : " << P->id << endl;
	P->faces.push_back(F->id);
	(F->voisins).push_back(P->id);
	if(F->BC == 1 || F->BC == -1)
	  (F->voisins).push_back(-1); //Car face au bord
	if(type == 5)
	  (F->c_reconstruction).push_back(0.5); //Maillage de Voronoi
	if(F->voisins.size() > 2) {
	  cout << "Numéro face : " << F->id << endl;
	  throw std::invalid_argument("Face a trop de voisins !");
	}
      }
    }

    /*if(F->voisins.size() == 1) {
      F->voisins.push_back(-1); // Cad face au bord
      F->BC = -1;
    }
    else
    F->BC = 0;*/

    int part_1 = F->voisins[0];
    int part_2 = -1;
    if(F->BC == 0)
      part_2 = F->voisins[1];

    //Vérification du sens de la normale
    if(part_1 != -1 && part_2 != -1) { //Face pas au bord
      if(Vector_3(solide[part_1].x0, solide[part_2].x0) * F->normale  < 0.)
	F->normale = -F->normale;
    }
    else if(part_1 == -1 || part_2 == -1) { //Face au bord
      Vector_3 bonne_direction = Vector_3(solide[F->voisins[0]].x0, F->centre);
      if(bonne_direction * F->normale < 0.)
	F->normale = -F->normale;
    }
  }


  //Calcul du tetrahèdre associé à chaque face pour le calcul du gradient
  /*if(type == 4) { //Seulement pour tetra
    for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){ //Boucle sur toutes les faces
      if(F->BC == 0) {
	//cout << "Face : " << F->id << endl;
	bool test = voisins_face(F->id);
	if(not(test)) {
	  //cout << "Face : " << F->id << " Pas de tetra associe a une face" << endl;
	  throw std::invalid_argument( "Pas de tetra associe a une face" );
	}
      }
    }
    }*/

  if(type == 4) { //Utilisation du Delaunay pour trouver les tétras associés à chaque face
    bool calcul_pre_traitement = false;
    std::ifstream pre_traitement("tetra_faces.txt"); //,ios::in);

    bool test = getline(pre_traitement, ligne);
    //cout << "test : " << test << endl;
    //Si pas de pretraitement
    if(test) {
      bool aux = true;
      for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
	if(aux && F->BC == 0) { //Importation des tetras associés aux faces internes
	  istringstream  stm(ligne);
	  int ele1,ele2,ele3,ele4;
	  double c1,c2,c3,c4;
	  stm >> ele1 >> ele2 >> ele3 >> ele4 >> c1 >> c2 >> c3 >> c4;
	  F->reconstruction.push_back(ele1);
	  F->reconstruction.push_back(ele2);
	  F->reconstruction.push_back(ele3);
	  F->reconstruction.push_back(ele4);
	  F->c_reconstruction.push_back(c1);
	  F->c_reconstruction.push_back(c2);
	  F->c_reconstruction.push_back(c3);
	  F->c_reconstruction.push_back(c4);
	  //cout << "face : " << F->id << endl;

	  aux = getline(pre_traitement, ligne);
	}
	/*else {
	  cout << "Fin fichier tetras associes faces !" << endl;
	  break;
	  }*/
      }
    }
    else {
      std::ifstream delaunay("delaunay.txt",ios::in);
      if(not(delaunay))
	throw std::invalid_argument( "Pas de terahedrisation de Delaunay !" );
      std::vector<Particule> tetra_delau; //Ensemble des particules contenant la tetra
      while(getline(delaunay, ligne)) { //Importation de la tetraedrisation de Delaunay
	istringstream  stm(ligne);
	int ele1,ele2,ele3,ele4;
	stm >> ele1 >> ele2 >> ele3 >> ele4;
	Particule P;
	P.vertices.push_back(ele1); //Numéros des Elements du solide qui forment chacun des tetra
	P.vertices.push_back(ele2);
	P.vertices.push_back(ele3); // - 1
	P.vertices.push_back(ele4);
	tetra_delau.push_back(P);
      }
      //cout << "ok Delaunay" << endl;
      std::ofstream tetra_faces("tetra_faces.txt",ios::out); //Sortie pour stocker tetras associés aux faces
      std::ofstream face_pb("face_pb.txt",ios::out); //Sorties pour les faces qui pose pb
      //Recherche du tetraèdre associé à chaque face
      for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
	bool tetra_ok = false;
	if(F->BC == 0) {
	  //cout << "Voisins : " << F->voisins[0] << " " << F->voisins[1] << endl;
	  for(std::vector<Particule>::iterator P=tetra_delau.begin();P!=tetra_delau.end();P++){
	    int part_1 = P->vertices[0];
	    int part_2 = P->vertices[1];
	    int voisin1 = P->vertices[2];
	    int voisin2 = P->vertices[3];
	  
	    double c1 = (Vector_3(solide[part_2].x0, F->centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_2].x0, solide[part_1].x0) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0) ));
	    double c2 = (Vector_3(solide[part_1].x0, F->centre) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_1].x0, solide[part_2].x0) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0) ));
	    double c3 = (Vector_3(solide[part_2].x0, F->centre) * cross_product(Vector_3(solide[part_2].x0, solide[part_1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_2].x0, solide[voisin1].x0) * cross_product(Vector_3(solide[part_2].x0, solide[part_1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0) ));
	    double c4 = (Vector_3(solide[part_2].x0, F->centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[part_1].x0)) ) / (Vector_3(solide[part_2].x0, solide[voisin2].x0) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[part_1].x0) ));

	    if( c1 >= 0. && c2 >= 0. && c3 >= 0. && c4 >= 0. && c1 < 1. && c2 < 1. && c3 < 1. && c4 < 1.) {
	      F->reconstruction.push_back(part_1);
	      F->reconstruction.push_back(part_2);
	      F->reconstruction.push_back(voisin1);
	      F->reconstruction.push_back(voisin2);
	      F->c_reconstruction.push_back(c1);
	      F->c_reconstruction.push_back(c2);
	      F->c_reconstruction.push_back(c3);
	      F->c_reconstruction.push_back(c4);
	      //cout << F->id << endl;
	      //if(c1 >= 0.&& c2 >= 0. && c3 >= 0.)
	      //cout << c1 << " " << c2 << " " << c3 << " " << c4 << " " << c1 + c2 + c3 + c4 - 1. << endl;
	      //cout << "face : " << F->id << " ok !" << endl;
	      tetra_ok = true;
	      //cout << "Face : " << F->id << " ok !" << endl;
	      break;
	    }
	  }
	  if( not(tetra_ok)) {
	    //cout << "Face : " << F->id;
	    //throw std::invalid_argument( " pas de tetra associe a une face !" );
	    /*cout << "Face : " << F->id << endl;
	      cout << "Voisins : " << F->voisins[0] << " " << F->voisins[1] << endl;
	      cout << "BC : " << F->BC << endl;
	      cout << "Num Vertex : " << F->vertex[0] << " " << F->vertex[1] << " " << F->vertex[2] << " " << endl;*/
	    face_pb << F->id << " " << vertex[F->vertex[0]].pos << " " << vertex[F->vertex[1]].pos << " " << vertex[F->vertex[2]].pos << endl;
	    //throw std::invalid_argument( " pas de tetra associe a la face !" );
	    bool test = voisins_face(F->id); //Dans ce cas, on va faire de l'extrapolation et utiliser l'ancienne méthode...
	    //cout << "face ok !" << endl;
	    if(not(test)) {
	      //cout << "Face : " << F->id << " Pas de tetra associe a une face" << endl;
	      throw std::invalid_argument( "Pas de tetra associ\'e a une face" );
	    }
	  }
	  /*else
	    cout << F->id << " : au bord" << endl;*/
      
	//Sortie contenant les éléments associés à chaque face pour faire le pré-traitement une seule fois
	//cout << "Enregistrement pret !" << endl;
	tetra_faces << F->reconstruction[0] << " " << F->reconstruction[1] << " " << F->reconstruction[2] << " " << F->reconstruction[3] << " " << F->c_reconstruction[0] << " " << F->c_reconstruction[1] << " " << F->c_reconstruction[2] << " " << F->c_reconstruction[3] << endl;
	}
      }
    }
  }
}

bool Solide::face_existe(Face f) { //Renvoie vraie si la face testée est déjà das faces
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
    if(*F == f) {
      return true;
    }
  }
  return false;
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
      double vol = std::abs(cross_product(Vector_3(solide[part_1].x0,solide[part_2].x0),Vector_3(solide[part_1].x0,solide[voisin1].x0))*Vector_3(solide[part_1].x0,solide[voisin2].x0)/6.); //Volume du tetra associé à la face
      if(vol > pow(10., -20.)) { //Attention avec cette valeur !! Il faut savoir qu'elle est là ! //-10.
	double c_part_1 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_2].x0, solide[part_1].x0) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0) ));
	double c_part_2 = (Vector_3(solide[part_1].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_1].x0, solide[part_2].x0) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0) ));
	double c_voisin1 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_2].x0, solide[part_1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_2].x0, solide[voisin1].x0) * cross_product(Vector_3(solide[part_2].x0, solide[part_1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0) ));
	double c_voisin2 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[part_1].x0)) ) / (Vector_3(solide[part_2].x0, solide[voisin2].x0) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[part_1].x0) ));

	//cout << "Calcul coords bary ok !" << endl;

	faces[num_face].c_reconstruction.push_back(c_part_1);
	faces[num_face].c_reconstruction.push_back(c_part_2);
	faces[num_face].c_reconstruction.push_back(c_voisin1);
	faces[num_face].c_reconstruction.push_back(c_voisin2);
        faces[num_face].reconstruction.push_back(part_1);
	faces[num_face].reconstruction.push_back(part_2);
	faces[num_face].reconstruction.push_back(voisin1);
        faces[num_face].reconstruction.push_back(voisin2);
	tetra_ok = true;
	break;	
      }
    }
  }
  return tetra_ok;
}

void Solide::Solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++)
    P->solve_position(dt, flag_2d, t, T);
}

void Solide::Solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++)
    P->solve_vitesse(dt, flag_2d, Amort, t , T);
    //P->solve_vitesse(dt, flag_2d, Amort, t , T, *this);
}

void Solide::Forces(const int& N_dim, const double& dt, const double& t, const double& T){
  Forces_internes(dt, t, T);
  /*for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++) {
    if(F->BC == -1)
      F->Fi = F->Fi; // + Forces_externes(t,T); //Forces ext s'appliquent sur les faces !
      }*/
}

void Solide::stresses(const double& t, const double& T){ //Calcul de la contrainte dans toutes les particules
  for(int i=0; i<faces.size(); i++){ //Calcul de la reconstruction sur chaque face
    if(faces[i].BC == 0) //cad face dans bulk et donc I_Dx reconstruit
      faces[i].I_Dx = Vector_3(0., 0., 0.); //Remise à zéro. Si particule sur le bord, on a bien I_Dx = (0., 0., 0.)
    //cout << "BC : " << faces[i].BC << endl;
    //Vector_3 test_pos(0., 0., 0.);
    if(faces[i].BC == 1) { //Dirichlet
      faces[i].I_Dx = Vector_3(0., 0., 0.);
      //faces[i].I_Dx = solide[faces[i].voisins[0]].Dx; //Dirichlet BC imposée fortement dans Mka ! old...
      //faces[i].I_Dx.vec[2] = faces[i].centre.z() /  3. * 4.;
      //cout << faces[i].I_Dx.vec[2] << endl;
      //faces[i].I_Dx = displacement_BC(faces[i].centre, faces[i].I_Dx, t, 0.); //Pour torsion
      //cout << "On impose bien les BC : " << faces[i].centre << " " << faces[i].I_Dx <<  endl;
      //if(faces[i].id == 0) //Test pour essayer de limiter la rotation...
	  //faces[i].I_Dx = Vector_3(0., 0., 0.);
      //double def_ref = 0.001 * t / T;
      //faces[i].I_Dx.vec[2] = faces[i].centre.z() * def_ref;
      //displacement_BC_bis(faces[i].centre, solide[faces[i].voisins[0]].Dx, t, 0.); //BC de Dirichlet
    }
    else if(faces[i].BC == -1) {
      faces[i].I_Dx = Vector_3(0., 0., 0.);
      //faces[i].I_Dx = solide[faces[i].voisins[0]].Dx; //Neumann
      //faces[i].I_Dx = displacement_BC(faces[i].centre, solide[faces[i].voisins[0]].Dx, t, 0.);
      //faces[i].I_Dx.vec[2] = displacement_BC_bis(faces[i].centre, solide[faces[i].voisins[0]].Dx, t, 0.);
    }
    else if(faces[i].BC == 0) { //Cad particule dans le bulk. Donc reconstruction !
      for(int j=0; j<faces[i].reconstruction.size() ; j++) {
	faces[i].I_Dx = faces[i].I_Dx + faces[i].c_reconstruction[j] * solide[faces[i].reconstruction[j]].Dx;
	//cout << "Interpolation : " << faces[i].I_Dx << endl;
      }
      //faces[i].I_Dx = (solide[faces[i].voisins[0]].Dx + solide[faces[i].voisins[1]].Dx) / 2.;

      /*if(abs(faces[i].I_Dx[2] - 4. / 3. * faces[i].centre.z()) > 0.00001)
	cout << "Problème reconstruction sur face : " << i << endl;*/
      /*if(sqrt((test_pos - Vector_3(Point_3(0.,0.,0.),faces[i].centre)).squared_length()) > pow(10., -11.))
	cout << "Problème reconstruction barycentre face : " << i << endl;*/
    }
  }
  
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    bool test_continuer;
    do {
      test_continuer = false;
      P->discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
      P->discrete_gradient.col2 = Vector_3(0., 0., 0.);
      P->discrete_gradient.col3 = Vector_3(0., 0., 0.);
      std::vector<int> num_faces; //Sert à récupérer le numéro des faces avec BC de Neumann si nécessaire
      for(int i=0 ; i < P->faces.size() ; i++){
        int f = P->faces[i];
	if(faces[f].BC != 0) { //car on ne fait ces calculs sur toutes les faces au bord
	  num_faces.push_back(f);
	}
	Vector_3 nIJ = faces[f].normale;
	if(faces[f].BC == 0 && nIJ * Vector_3(P->x0, faces[f].centre) < 0.)
	  nIJ = -nIJ; //Normale pas dans le bon sens...
	Matrix Dij_n(tens_sym(faces[f].I_Dx,  nIJ));
	if(faces[f].BC == 0) //On ajoute pas les faces de Neumann ni de Dirichlet car on doit recalculer la valeur sur la face// >= 0 avant test
	  P->discrete_gradient += faces[f].S /  P->V * Dij_n;
	//else if(faces[f].BC == 1) //On ajoute la partie Dirichlet imposée. Si toute la face était en dirichlet elle finirait comme une face intérieure au-dessus...
	  //P->discrete_gradient += -faces[f].S /  P->V * ((Dij_n * faces[F].normale) * faces[F].normale) * faces[F].normale;
      }

      P->contrainte = lambda * (P->discrete_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_gradient - P->epsilon_p); //Premier calcul pour calcul des déplacements sur bords de Neumann
      
      //On reconstruit la valeur du déplacement sur les faces de Neumann
      if(num_faces.size() > 0) {
	reconstruction_faces_neumann(num_faces, P->contrainte, t, P->V, T); //Calcul la valeur des déplacements sur faces de Neumann
	
	//On recalcul le gradient discret
	P->discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
	P->discrete_gradient.col2 = Vector_3(0., 0., 0.);
	P->discrete_gradient.col3 = Vector_3(0., 0., 0.);
	for(int i=0 ; i < P->faces.size() ; i++){ //On recalcule les contraintes avec toutes les contributions
	  int f = P->faces[i];

	  //On se fout de la reconstruction, on a des déplacements imposés !
	  /*double def_ref = 3. / 4. * t / T;
	  faces[f].I_Dx.vec[2] = faces[f].centre.z() * def_ref;
	  faces[f].I_Dx.vec[0] = -0.3 * faces[f].centre.x() * def_ref;
	  faces[f].I_Dx.vec[1] = -0.3 * faces[f].centre.y() * def_ref; //On impose les positions pour le test
	  */
	  
	  Vector_3 nIJ = faces[f].normale;
	  if(faces[f].BC == 0 && nIJ * Vector_3(P->x0, faces[f].centre) < 0.)
	    nIJ = -nIJ; //Normale pas dans le bon sens...
	  Matrix Dij_n(tens_sym(faces[f].I_Dx,  nIJ));
	  P->discrete_gradient += faces[f].S /  P->V * Dij_n;
	}
      }
      
      P->contrainte = lambda * (P->discrete_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_gradient - P->epsilon_p); //Premier calcul pour calcul des déplacements sur bords de Neumann
      /*for(int i=0 ; i < P->faces.size() ; i++){ //On recalcule les contraintes avec toutes les contributions
	  int f = P->faces[i];
	  if(faces[f].BC == -1 && sqrt((P->contrainte * faces[f].normale).squared_length()) > 1.)
	    cout << "Problème Neumann homogène : " << sqrt((P->contrainte * faces[f].normale).squared_length()) << ". Num face : " << faces[f].id << endl;
	    }*/

      //Plastification si on dépasse le critère  
      P->seuil_elas = A; // + B * pow(P->def_plas_cumulee, n);
      if((P->contrainte - H * P->epsilon_p).VM() > P->seuil_elas) { //On sort du domaine élastique.
	//while( (P->contrainte).VM() > A) { //On dépasse le critère plastique on refait un return mapping pour essayer de converger
	//Plastification
	Matrix n_elas( 1. / ((P->contrainte).dev()).norme() * (P->contrainte).dev() ); //Normale au domaine élastique de Von Mises
	double delta_p = ((P->contrainte - H * P->epsilon_p).VM() - A) / (2*mu + H);
	P->def_plas_cumulee += delta_p;
	P->epsilon_p += delta_p * n_elas;
	P->contrainte = lambda * (P->discrete_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_gradient - P->epsilon_p); //Recalcul des contraintes après plastification
	}
      if(num_faces.size() > 0) { //On vérifie qu'on a toujours les bonnes BC de Neumann
	for(int i=0 ; i < num_faces.size() ; i++) {
	  if( sqrt( (P->contrainte * faces[i].normale).squared_length()) > 1.) {
	    test_continuer = true;
	    break;
	  }
	}
      }
    } while((P->contrainte - H * P->epsilon_p).VM() > P->seuil_elas && test_continuer); //Ajouter dans le test que les conditions de bord doivent être respectées ?
    /*if((P->contrainte).VM() > P->seuil_elas)
      cout << "Von Mises : " << (P->contrainte - H * P->epsilon_p).VM() << endl;*/
  }
}

void Solide::reconstruction_faces_neumann(std::vector<int> num_faces, const Matrix& contrainte, const double& t, const double& V, const double& T) {
  if(num_faces.size() == 1){  //Solution directe sans inversion de matrice
    //Reconstruction de la valeur sur face de Neumann Homogène s'il y en a une
    int F = num_faces[0];
    if(faces[F].BC == -1) { //Cas Neumann Complet ==-1
      //Test avec inversion de la matrice
      Eigen::Matrix<double, 3, 1> b; //Vecteur second membre. Neumann homogène pour l'instant
      Eigen::Matrix<double, 3, 1> x; //Contient les valeurs aux faces
      Eigen::MatrixXd Mat(3,3); //Premier bloc diagonal
      
      Mat << (lambda + mu) * faces[F].normale.x() * faces[F].normale.x() + mu, (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.y() * faces[F].normale.y() + mu, (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.z() * faces[F].normale.z() + mu;
      Mat *= faces[F].S / V;

      b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.);
      //b *= 2.; //Semble marcher. Pk ????

      //Inversion du système !
      typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
      Eigen::FullPivLU<Matrix3x3> lu(Mat);
      if( lu.rank() == 3) //Test voir si système inversible...
	x = Mat.lu().solve(b); //Problème avec les valeurs de x !!!!
      else { //Calcul de la pseudo-inverse pour minimisation de l'écart aux moindres carrés.
	Eigen::CompleteOrthogonalDecomposition<Matrix3x3> mat(Mat);
	x = mat.solve(b);
      }

      /*double def_ref = 0.001 * t / T;
      cout << "Attendu : " << faces[F].centre.z() * def_ref << " " << -0.3 * faces[F].centre.x() * def_ref << " " << -0.3 * faces[F].centre.y() * def_ref << endl;
      cout << "Deplacement normal : " << x(2) << endl;
      cout << "Deplacements tangents : " << x(0) << " " << x(1) << endl;*/
      //double cc1 = 2. * x(0); double cc2 = 2. * x(1); double cc3 = 2. * x(2);
      
      faces[F].I_Dx.vec[0] = x(0); faces[F].I_Dx.vec[1] = x(1); faces[F].I_Dx.vec[2] = x(2);
      //faces[F].I_Dx.vec[0] = cc1; faces[F].I_Dx.vec[1] = cc2; faces[F].I_Dx.vec[2] = cc3;

      //cout << "Valeur interpolée : "  << faces[F].I_Dx << endl;
    }
    else if(faces[F].BC == 1) { //Calcul direct avec vecteurs tangents
      /*Eigen::Matrix<double, 3, 1> b; //Vecteur second membre. Neumann homogène pour l'instant
      Eigen::Matrix<double, 3, 1> x; //Contient les valeurs aux faces
      Eigen::MatrixXd Mat(3,3); //Premier bloc diagonal
      
      Mat << (lambda + mu) * faces[F].normale.x() * faces[F].normale.x() + mu, (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.y() * faces[F].normale.y() + mu, (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.z() * faces[F].normale.z() + mu;
      Mat *= faces[F].S / V;
      Mat(2,0) = 0.;
      Mat(2,1) = 0.;
      Mat(2,2) = 1.;*/

      Eigen::Matrix<double, 4, 1> b; //4 //Vecteur second membre. Neumann homogène pour l'instant
      Eigen::Matrix<double, 3, 1> x; //3 //Contient les valeurs aux faces
      Eigen::MatrixXd Mat(4,3); //4 //Premier bloc diagonal
				    
      Mat << (lambda + mu) * faces[F].normale.x() * faces[F].normale.x() + mu, (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.y() * faces[F].normale.y() + mu, (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.z() * faces[F].normale.z() + mu, 0., 0., 0.;
      Mat *= faces[F].S / V;
      double nrm = Mat.norm();
      Mat(3,2) = nrm; //1.e17 //Attention la valeur mise ici doit coller avec celle dans second membre et ordre de grandeur des autres valeurs dans la matrice !
      //Mat(3,2) = 1.; //Attention la valeur mise ici doit coller avec celle dans second membre et ordre de grandeur des autres valeurs dans la matrice !


      //cout << "Matrice inversée : " << Mat << endl;

      double def_ref = 0.001 * t / T;
      //b << (-contrainte * faces[F].normale + ((contrainte * faces[F].normale) * faces[F].normale) * faces[F].normale + faces[F].S / V * (lambda + 2.*mu) * faces[F].centre.z() * def_ref * faces[F].normale) * Vector_3(1.,0.,0.), (-contrainte * faces[F].normale + ((contrainte * faces[F].normale) * faces[F].normale) * faces[F].normale  + faces[F].S / V * (lambda + 2.*mu) * faces[F].centre.z() * def_ref * faces[F].normale) * Vector_3(0.,1.,0.), (-contrainte * faces[F].normale + ((contrainte * faces[F].normale) * faces[F].normale) * faces[F].normale  + faces[F].S / V * (lambda + 2.*mu) * faces[F].centre.z() * def_ref * faces[F].normale) * Vector_3(0.,0.,1.), /*nrm * */faces[F].centre.z() * def_ref; //displacement_BC_bis(faces[F].centre, solide[faces[F].voisins[0]].Dx, t, 0.); //BC de Dirichlet

      b << (-contrainte * faces[F].normale) * Vector_3(1.,0.,0.), (-contrainte * faces[F].normale) * Vector_3(0.,1.,0.), (-contrainte * faces[F].normale) * Vector_3(0.,0.,1.), nrm * faces[F].centre.z() * def_ref; //displacement_BC_bis(faces[F].centre, solide[faces[F].voisins[0]].Dx, t, 0.); //

      //cout << "Second membre : " << b << endl;

      //Inversion du système !
      /*typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
      Eigen::FullPivLU<Matrix3x3> lu(Mat);
      if( lu.rank() == 3) //Test voir si système inversible...
	x = Mat.lu().solve(b); //Problème avec les valeurs de x !!!!
	else {*/
	//Calcul de la pseudo-inverse pour minimisation de l'écart aux moindres carrés.
      typedef Eigen::Matrix<double, 4, 3> Matrix4x3;
      Eigen::CompleteOrthogonalDecomposition<Matrix4x3> mat(Mat);
      x = mat.solve(b);
	//}

      //cout << "Attendu : " << faces[F].centre.z() * def_ref << endl; //displacement_BC_bis(faces[F].centre, solide[faces[F].voisins[0]].Dx, t, 0.) << endl; // << " " << -0.3 * faces[F].centre.x() * def_ref << " " << -0.3 * faces[F].centre.y() * def_ref << endl;
      //cout << "Deplacement normal : " << x(2) << endl;
      //cout << "Deplacements tangents : " << x(0) << " " << x(1) << endl;
      
      faces[F].I_Dx.vec[0] = x(0); faces[F].I_Dx.vec[1] = x(1); faces[F].I_Dx.vec[2] = x(2);
    }

    //Test contraintes
    /*Vector_3 nIJ = faces[F].normale;
    Matrix Dij_n(tens_sym(faces[F].I_Dx,  nIJ)); //- P->Dx
    P->discrete_gradient += faces[F].S /  P->V * Dij_n;
    
    contrainte = lambda * (P->discrete_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_gradient - P->epsilon_p); //Calcul des contraintes complètes
    if(sqrt((contrainte * faces[F].normale).squared_length()) > 0.0001)
    cout << "Contrainte sur bord 1 de Neumann : " << sqrt((contrainte * faces[F].normale).squared_length()) << endl;*/
  }
  else if(num_faces.size() == 2) { //Inversion d'un système linéaire de 6 équations avec Eigen
    //cout << "2 faces sur bord de Neumann !" << endl;
    int F = num_faces[0];
    int Fp = num_faces[1];
    Eigen::Matrix<double, 6, 1> x; //Contient les valeurs aux faces

    Eigen::MatrixXd A_FF(3,3); //Premier bloc diagonal
    A_FF << (lambda + mu) * faces[F].normale.x() * faces[F].normale.x() + mu, (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.y() * faces[F].normale.y() + mu, (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.z() * faces[F].normale.z() + mu;
    A_FF *= faces[F].S / V;

    Eigen::MatrixXd A_FFp(3,3); //Premier bloc hors-diagonale
    A_FFp << (lambda + mu) * faces[Fp].normale.x() * faces[F].normale.x() + mu * (faces[F].normale * faces[Fp].normale), lambda * faces[Fp].normale.y() * faces[F].normale.x() + mu * faces[Fp].normale.x() * faces[F].normale.y(),   lambda * faces[Fp].normale.z() * faces[F].normale.x() + mu * faces[Fp].normale.x() * faces[F].normale.z(), lambda * faces[Fp].normale.x() * faces[F].normale.y() + mu * faces[F].normale.x() * faces[Fp].normale.y(),  (lambda + mu) * faces[Fp].normale.y() * faces[F].normale.y() + mu * (faces[F].normale * faces[Fp].normale), lambda * faces[F].normale.y() * faces[Fp].normale.z() + mu * faces[Fp].normale.y() * faces[F].normale.z(),  lambda * faces[Fp].normale.x() * faces[F].normale.z() + mu * faces[F].normale.x() * faces[Fp].normale.z(),  lambda * faces[Fp].normale.y() * faces[F].normale.z() + mu * faces[F].normale.y() * faces[Fp].normale.z(), (lambda + mu) * faces[Fp].normale.z() * faces[F].normale.z() + mu * (faces[F].normale * faces[Fp].normale);
    A_FFp *= faces[Fp].S / V;

    Eigen::MatrixXd A_FpF(3,3); //Second bloc hors-diagonale
    A_FpF << (lambda + mu) * faces[F].normale.x() * faces[Fp].normale.x() + mu * (faces[F].normale * faces[Fp].normale), lambda * faces[Fp].normale.x() * faces[F].normale.y() + mu * faces[F].normale.x() * faces[Fp].normale.y(),   lambda * faces[Fp].normale.x() * faces[F].normale.z() + mu * faces[F].normale.x() * faces[Fp].normale.z(), lambda * faces[F].normale.x() * faces[Fp].normale.y() + mu * faces[Fp].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[Fp].normale.y() * faces[F].normale.y() + mu * (faces[F].normale * faces[Fp].normale), lambda * faces[Fp].normale.y() * faces[F].normale.z() + mu * faces[F].normale.y() * faces[Fp].normale.z(),  lambda * faces[F].normale.x() * faces[Fp].normale.z() + mu * faces[Fp].normale.x() * faces[F].normale.z(),  lambda * faces[F].normale.y() * faces[Fp].normale.z() + mu * faces[Fp].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[Fp].normale.z() * faces[F].normale.z() + mu * (faces[F].normale * faces[Fp].normale);
    A_FpF *= faces[F].S / V;

    Eigen::MatrixXd A_FpFp(3,3); //Second bloc diagonal
    A_FpFp << (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.x() + mu, (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.y(),  (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.z(),  (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.y(),  (lambda + mu) * faces[Fp].normale.y() * faces[Fp].normale.y() + mu, (lambda + mu) * faces[Fp].normale.y() * faces[Fp].normale.z(), (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.z(), (lambda + mu) * faces[Fp].normale.y() * faces[Fp].normale.z(), (lambda + mu) * faces[Fp].normale.z() * faces[Fp].normale.z() + mu;
    A_FpFp *= faces[Fp].S / V;

    if(faces[F].BC == -1 && faces[Fp].BC == -1) {
      Eigen::Matrix<double, 6, 6> Mat; //Matrice à inverser
      Eigen::Matrix<double, 6, 1> b; //Vecteur second membre. Neumann homogène pour l'instant

      //Assemblage de la matrice
      Mat.topLeftCorner<3,3>() = A_FF;
      Mat.topRightCorner<3,3>() = A_FFp;
      Mat.bottomLeftCorner<3,3>() = A_FpF;
      Mat.bottomRightCorner<3,3>() = A_FpFp;

    //Assemblage du second membre
    b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.),  ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,0.,1.);

    typedef Eigen::Matrix<double, 6, 6> Matrix6x6;
    Eigen::FullPivLU<Matrix6x6> lu(Mat);
    if( lu.rank() == 6) //Test voir si système inversible...
      x = Mat.lu().solve(b); //Problème avec les valeurs de x !!!!
    else { //Calcul de la pseudo-inverse pour minimisation de l'écart aux moindres carrés.
      typedef Eigen::Matrix<double, 6, 6> Matrix6x6;
      Eigen::CompleteOrthogonalDecomposition<Matrix6x6> mat(Mat);
      x = mat.solve(b);
      }

    faces[F].I_Dx.vec[0] = x(0); faces[F].I_Dx.vec[1] = x(1); faces[F].I_Dx.vec[2] = x(2); //Première face de Neumann
    faces[Fp].I_Dx.vec[0] = x(3); faces[Fp].I_Dx.vec[1] = x(4); faces[Fp].I_Dx.vec[2] = x(5); //Deuxième face de Neumann
    }
    else if(faces[F].BC == 1) {
      Eigen::Matrix<double, 7, 6> Mat; //Matrice à inverser
      Eigen::Matrix<double, 7, 1> b; //Vecteur second membre. Neumann homogène pour l'instant

      //Assemblage de la matrice
      /*Mat.topLeftCorner<3,3>() = A_FF;
      Mat.topRightCorner<3,3>() = A_FFp;
      Mat.bottomLeftCorner<3,3>() = A_FpF;
      Mat.bottomRightCorner<3,3>() = A_FpFp;*/
      Mat.block<3,3>(0,0) = A_FF;
      Mat.block<3,3>(0,3) = A_FFp;
      Mat.block<3,3>(4,0) = A_FpF;
      Mat.block<3,3>(4,3) = A_FpFp;

      //Assemblage du second membre
      //b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.),  ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,0.,1.);
      Mat(3,0) = 0.;
      Mat(3,1) = 0.;
      //Mat(3,2) = Mat(0,0); //Attention la valeur mise ici doit coller avec celle dans second membre et ordre de grandeur des autres valeurs dans la matrice !
      double nrm = Mat.norm();
      Mat(3,2) = nrm; //1. //1.e15 //Attention la valeur mise ici doit coller avec celle dans second membre et ordre de grandeur des autres valeurs dans la matrice !
      Mat(3,3) = 0.;
      Mat(3,4) = 0.;
      Mat(3,5) = 0.;
      //cout << "Matrice à inverser : " << Mat << endl;
      double def_ref = 0.001 * t / T;
      //b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), displacement_BC_bis(faces[F].centre, solide[faces[F].voisins[0]].Dx, t, 0.), ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,0.,1.);
      b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.), nrm * faces[F].centre.z() * def_ref/* displacement_BC_bis(faces[F].centre, solide[faces[F].voisins[0]].Dx, t, 0.)*/, ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,0.,1.);

      //Résolution
      typedef Eigen::Matrix<double, 7, 6> Matrix7x6;
      Eigen::CompleteOrthogonalDecomposition<Matrix7x6> mat(Mat);
      x = mat.solve(b);

      //cout << "Attendu : " << faces[F].centre.z() * def_ref << endl;//displacement_BC_bis(faces[F].centre, solide[faces[F].voisins[0]].Dx, t, 0.) << endl; //faces[F].centre.z() * def_ref << " " << -0.3 * faces[F].centre.x() * def_ref << " " << -0.3 * faces[F].centre.y() * def_ref << endl;
      //cout << "Deplacement normal : " << x(2) << endl;
      //cout << "Attendu : " << -0.3 * faces[F].centre.x() * def_ref << " " << -0.3 * faces[F].centre.y() * def_ref << endl;
      //cout << "Deplacements bord : " << x(0) << " " << x(1) << " et " << x(3) << " " << x(4) << endl;

      faces[F].I_Dx.vec[0] = x(0); faces[F].I_Dx.vec[1] = x(1); faces[F].I_Dx.vec[2] = x(2); //Première face de Neumann
      faces[Fp].I_Dx.vec[0] = x(3); faces[Fp].I_Dx.vec[1] = x(4); faces[Fp].I_Dx.vec[2] = x(5); //Deuxième face de Neumann
    }
    else if(faces[Fp].BC == 1) {
      Eigen::Matrix<double, 7, 6> Mat; //Matrice à inverser
      Eigen::Matrix<double, 7, 1> b; //Vecteur second membre. Neumann homogène pour l'instant

      //Assemblage de la matrice
      /*Mat.topLeftCorner<3,3>() = A_FF;
      Mat.topRightCorner<3,3>() = A_FFp;
      Mat.bottomLeftCorner<3,3>() = A_FpF;
      Mat.bottomRightCorner<3,3>() = A_FpFp;*/
      Mat.block<3,3>(0,0) = A_FF;
      Mat.block<3,3>(0,3) = A_FFp;
      Mat.block<3,3>(4,0) = A_FpF;
      Mat.block<3,3>(4,3) = A_FpFp;
      Mat(3,0) = 0.;
      Mat(3,1) = 0.;
      Mat(3,2) = 0.;
      Mat(3,3) = 0.;
      Mat(3,4) = 0.;
      //cout << "Matrice à inverser : " << Mat << endl;
      //Mat(3,5) = Mat(0,0); //Attention la valeur mise ici doit coller avec celle dans second membre et ordre de grandeur des autres valeurs dans la matrice !
      double nrm = Mat.norm();
      Mat(3,5) = nrm; //1. //1.e15 //Attention la valeur mise ici doit coller avec celle dans second membre et ordre de grandeur des autres valeurs dans la matrice !
      double def_ref = 0.001 * t / T;
      //b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.),  ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.),  displacement_BC_bis(faces[Fp].centre, solide[faces[Fp].voisins[0]].Dx, t, 0.);
      b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.),  nrm * faces[F].centre.z() * def_ref /*displacement_BC_bis(faces[Fp].centre, solide[faces[Fp].voisins[0]].Dx, t, 0.)*/, ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,0.,1.);

      //Résolution
      typedef Eigen::Matrix<double, 7, 6> Matrix7x6;
      Eigen::CompleteOrthogonalDecomposition<Matrix7x6> mat(Mat);
      x = mat.solve(b);
      /*cout << "Attendu : " << -0.3 * faces[F].centre.x() * def_ref << " " << -0.3 * faces[F].centre.y() * def_ref << endl;
      //cout << "Attendu : " << faces[F].centre.z() * def_ref << endl;
      cout << "Deplacements bord : " << x(0) << " " << x(1) << " et " << x(3) << " " << x(4) << endl;*/

      //cout << "Attendu : " << faces[F].centre.z() * def_ref << endl; //displacement_BC_bis(faces[Fp].centre, solide[faces[Fp].voisins[0]].Dx, t, 0.) << endl; //faces[F].centre.z() * def_ref << " " << -0.3 * faces[F].centre.x() * def_ref << " " << -0.3 * faces[F].centre.y() * def_ref << endl;
      //cout << "Deplacement normal : " << x(5) << endl;
      //cout << "Attendu : " << -0.3 * faces[F].centre.x() * def_ref << " " << -0.3 * faces[F].centre.y() * def_ref << endl;
      //cout << "Deplacements bord : " << x(0) << " " << x(1) << " et " << x(3) << " " << x(4) << endl;

      faces[F].I_Dx.vec[0] = x(0); faces[F].I_Dx.vec[1] = x(1); faces[F].I_Dx.vec[2] = x(2); //Première face de Neumann
      faces[Fp].I_Dx.vec[0] = x(3); faces[Fp].I_Dx.vec[1] = x(4); faces[Fp].I_Dx.vec[2] = x(5); //Deuxième face de Neumann
    }
      
    //Inversion du système !
    /*typedef Eigen::Matrix<double, 6, 6> Matrix6x6;
    Eigen::FullPivLU<Matrix6x6> lu(Mat);
    if( lu.rank() == 6) //Test voir si système inversible...
      x = Mat.lu().solve(b); //Problème avec les valeurs de x !!!!
    else { //Calcul de la pseudo-inverse pour minimisation de l'écart aux moindres carrés.
      //cout << "Rang : " << lu.rank() << endl;
      //cout << faces[F].id << " " << faces[Fp].id << " : Neumann pas inversible !" << endl;
      typedef Eigen::Matrix<double, 7, 6> Matrix7x6;
      Eigen::CompleteOrthogonalDecomposition<Matrix7x6> mat(Mat);
      x = mat.solve(b);
      }*/

    //faces[F].I_Dx.vec[0] = x(0); faces[F].I_Dx.vec[1] = x(1); faces[F].I_Dx.vec[2] = x(2); //Première face de Neumann
    //faces[Fp].I_Dx.vec[0] = x(3); faces[Fp].I_Dx.vec[1] = x(4); faces[Fp].I_Dx.vec[2] = x(5); //Deuxième face de Neumann
    //cout << "Deplacements bord : " << x(3) << " et " << x(5) << endl;

    //Test pour voir si c'est bon ici ou pas...
    /*contrainte = lambda * (P->discrete_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_gradient - P->epsilon_p); //Calcul des contraintes complètes
    if(sqrt((contrainte * faces[F].normale).squared_length()) > 0.0001)
      cout << "Contrainte sur bord 1 de Neumann : " << sqrt((contrainte * faces[F].normale).squared_length()) << endl;
      if(sqrt((contrainte * faces[Fp].normale).squared_length()) > 0.0001)
      cout << "Contrainte sur bord 2 de Neumann : " << sqrt((contrainte * faces[Fp].normale).squared_length()) << endl;*/

  }
  else if(num_faces.size() == 3) { //Inversion d'un système linéaire de 9 équations avec Eigen
    //Il va falloir coder la partie pour conditions mixtes Neumann-Dirichlet !
    cout << "Particule à 3 faces de Neumann !" << endl; //Reprendre l'écriture des matrices !
    int F = num_faces[0];
    int Fp = num_faces[1];
    int Fpp = num_faces[2];
    Eigen::Matrix<double, 9, 9> Mat; //Matrice à inverser
    Eigen::Matrix<double, 9, 1> b; //Vecteur second membre. Neumann homogène pour l'instant
    Eigen::Matrix<double, 9, 1> x; //Contient les valeurs aux faces

    Eigen::MatrixXd A_FF(3,3); //Premier bloc diagonal
    A_FF << (lambda + mu) * faces[F].normale.x() * faces[F].normale.x() + mu, (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(),  (lambda + mu) * faces[F].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[F].normale.y() * faces[F].normale.y() + mu, (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.x() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[F].normale.z() * faces[F].normale.z() + mu;
    A_FF *= faces[F].S / V;

    Eigen::MatrixXd A_FFp(3,3); //Premier bloc hors-diagonale
    A_FFp << (lambda + mu) * faces[Fp].normale.x() * faces[F].normale.x() + mu * (faces[F].normale * faces[Fp].normale), lambda * faces[Fp].normale.y() * faces[F].normale.x() + mu * faces[Fp].normale.x() * faces[F].normale.y(),   lambda * faces[Fp].normale.z() * faces[F].normale.x() + mu * faces[Fp].normale.x() * faces[F].normale.z(), lambda * faces[Fp].normale.x() * faces[F].normale.y() + mu * faces[F].normale.x() * faces[Fp].normale.y(),  (lambda + mu) * faces[Fp].normale.y() * faces[F].normale.y() + mu * (faces[F].normale * faces[Fp].normale), lambda * faces[F].normale.y() * faces[Fp].normale.z() + mu * faces[Fp].normale.y() * faces[F].normale.z(),  lambda * faces[Fp].normale.x() * faces[F].normale.z() + mu * faces[F].normale.x() * faces[Fp].normale.z(),  lambda * faces[Fp].normale.y() * faces[F].normale.z() + mu * faces[F].normale.y() * faces[Fp].normale.z(), (lambda + mu) * faces[Fp].normale.z() * faces[F].normale.z() + mu * (faces[F].normale * faces[Fp].normale);
    A_FFp *= faces[Fp].S / V;

    Eigen::MatrixXd A_FFpp(3,3); //Second bloc hors-diagonale
    A_FFpp << (lambda + mu) * faces[Fpp].normale.x() * faces[F].normale.x() + mu * (faces[F].normale * faces[Fpp].normale), lambda * faces[Fpp].normale.y() * faces[F].normale.x() + mu * faces[Fpp].normale.x() * faces[F].normale.y(),   lambda * faces[Fpp].normale.z() * faces[F].normale.x() + mu * faces[Fpp].normale.x() * faces[F].normale.z(), lambda * faces[Fpp].normale.x() * faces[F].normale.y() + mu * faces[F].normale.x() * faces[Fpp].normale.y(),  (lambda + mu) * faces[Fpp].normale.y() * faces[F].normale.y() + mu * (faces[F].normale * faces[Fpp].normale), lambda * faces[F].normale.y() * faces[Fpp].normale.z() + mu * faces[Fpp].normale.y() * faces[F].normale.z(),  lambda * faces[Fpp].normale.x() * faces[F].normale.z() + mu * faces[F].normale.x() * faces[Fpp].normale.z(),  lambda * faces[Fpp].normale.y() * faces[F].normale.z() + mu * faces[F].normale.y() * faces[Fpp].normale.z(), (lambda + mu) * faces[Fpp].normale.z() * faces[F].normale.z() + mu * (faces[F].normale * faces[Fpp].normale);
    A_FFpp *= faces[Fpp].S / V;

    Eigen::MatrixXd A_FpF(3,3); //Troisième bloc hors-diagonale
    A_FpF << (lambda + mu) * faces[F].normale.x() * faces[Fp].normale.x() + mu * (faces[F].normale * faces[Fp].normale), lambda * faces[Fp].normale.x() * faces[F].normale.y() + mu * faces[F].normale.x() * faces[Fp].normale.y(),   lambda * faces[Fp].normale.x() * faces[F].normale.z() + mu * faces[F].normale.x() * faces[Fp].normale.z(), lambda * faces[F].normale.x() * faces[Fp].normale.y() + mu * faces[Fp].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[Fp].normale.y() * faces[F].normale.y() + mu * (faces[F].normale * faces[Fp].normale), lambda * faces[Fp].normale.y() * faces[F].normale.z() + mu * faces[F].normale.y() * faces[Fp].normale.z(),  lambda * faces[F].normale.x() * faces[Fp].normale.z() + mu * faces[Fp].normale.x() * faces[F].normale.z(),  lambda * faces[F].normale.y() * faces[Fp].normale.z() + mu * faces[Fp].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[Fp].normale.z() * faces[F].normale.z() + mu * (faces[F].normale * faces[Fp].normale);
    A_FpF *= faces[F].S / V;

    Eigen::MatrixXd A_FpFp(3,3); //Second bloc diagonal
    A_FpFp << (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.x() + mu, (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.y(),  (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.z(),  (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.y(),  (lambda + mu) * faces[Fp].normale.y() * faces[Fp].normale.y() + mu, (lambda + mu) * faces[Fp].normale.y() * faces[Fp].normale.z(), (lambda + mu) * faces[Fp].normale.x() * faces[Fp].normale.z(), (lambda + mu) * faces[Fp].normale.y() * faces[Fp].normale.z(), (lambda + mu) * faces[Fp].normale.z() * faces[Fp].normale.z() + mu;
    A_FpFp *= faces[Fp].S / V;

    Eigen::MatrixXd A_FpFpp(3,3); //Quatrième bloc hors-diagonale
    A_FpFpp << (lambda + mu) * faces[Fpp].normale.x() * faces[Fp].normale.x() + mu * (faces[Fpp].normale * faces[Fp].normale), lambda * faces[Fp].normale.x() * faces[Fpp].normale.y() + mu * faces[Fpp].normale.x() * faces[Fp].normale.y(),   lambda * faces[Fp].normale.x() * faces[Fpp].normale.z() + mu * faces[Fpp].normale.x() * faces[Fp].normale.z(), lambda * faces[Fpp].normale.x() * faces[Fp].normale.y() + mu * faces[Fp].normale.x() * faces[Fpp].normale.y(),  (lambda + mu) * faces[Fp].normale.y() * faces[Fpp].normale.y() + mu * (faces[Fpp].normale * faces[Fp].normale), lambda * faces[Fp].normale.y() * faces[Fpp].normale.z() + mu * faces[Fpp].normale.y() * faces[Fp].normale.z(),  lambda * faces[Fpp].normale.x() * faces[Fp].normale.z() + mu * faces[Fp].normale.x() * faces[Fpp].normale.z(),  lambda * faces[Fpp].normale.y() * faces[Fp].normale.z() + mu * faces[Fp].normale.y() * faces[Fpp].normale.z(), (lambda + mu) * faces[Fp].normale.z() * faces[Fpp].normale.z() + mu * (faces[Fpp].normale * faces[Fp].normale);
    A_FpFpp *= faces[Fpp].S / V;

    Eigen::MatrixXd A_FppF(3,3); //Cinquième bloc hors-diagonale
    A_FppF << (lambda + mu) * faces[F].normale.x() * faces[Fpp].normale.x() + mu * (faces[F].normale * faces[Fpp].normale), lambda * faces[Fpp].normale.x() * faces[F].normale.y() + mu * faces[F].normale.x() * faces[Fpp].normale.y(),   lambda * faces[Fpp].normale.x() * faces[F].normale.z() + mu * faces[F].normale.x() * faces[Fpp].normale.z(), lambda * faces[F].normale.x() * faces[Fpp].normale.y() + mu * faces[Fpp].normale.x() * faces[F].normale.y(),  (lambda + mu) * faces[Fpp].normale.y() * faces[F].normale.y() + mu * (faces[F].normale * faces[Fpp].normale), lambda * faces[Fpp].normale.y() * faces[F].normale.z() + mu * faces[F].normale.y() * faces[Fpp].normale.z(),  lambda * faces[F].normale.x() * faces[Fpp].normale.z() + mu * faces[Fpp].normale.x() * faces[F].normale.z(),  lambda * faces[F].normale.y() * faces[Fpp].normale.z() + mu * faces[Fpp].normale.y() * faces[F].normale.z(), (lambda + mu) * faces[Fpp].normale.z() * faces[F].normale.z() + mu * (faces[F].normale * faces[Fpp].normale);
    A_FppF *= faces[F].S / V;

    Eigen::MatrixXd A_FppFp(3,3); //Sixième bloc hors-diagonale
    A_FppFp << (lambda + mu) * faces[Fp].normale.x() * faces[Fpp].normale.x() + mu * (faces[Fp].normale * faces[Fpp].normale), lambda * faces[Fpp].normale.x() * faces[Fp].normale.y() + mu * faces[Fp].normale.x() * faces[Fpp].normale.y(),   lambda * faces[Fpp].normale.x() * faces[Fp].normale.z() + mu * faces[Fp].normale.x() * faces[Fpp].normale.z(), lambda * faces[Fp].normale.x() * faces[Fpp].normale.y() + mu * faces[Fpp].normale.x() * faces[Fp].normale.y(),  (lambda + mu) * faces[Fpp].normale.y() * faces[Fp].normale.y() + mu * (faces[Fp].normale * faces[Fpp].normale), lambda * faces[Fpp].normale.y() * faces[Fp].normale.z() + mu * faces[Fp].normale.y() * faces[Fpp].normale.z(),  lambda * faces[Fp].normale.x() * faces[Fpp].normale.z() + mu * faces[Fpp].normale.x() * faces[Fp].normale.z(),  lambda * faces[Fp].normale.y() * faces[Fpp].normale.z() + mu * faces[Fpp].normale.y() * faces[Fp].normale.z(), (lambda + mu) * faces[Fpp].normale.z() * faces[Fp].normale.z() + mu * (faces[Fp].normale * faces[Fpp].normale);
    A_FppFp *= faces[Fp].S / V;

    Eigen::MatrixXd A_FppFpp(3,3); //Troisième bloc diagonal
    A_FppFpp << (lambda + mu) * faces[Fpp].normale.x() * faces[Fpp].normale.x() + mu, (lambda + mu) * faces[Fpp].normale.x() * faces[Fpp].normale.y(),  (lambda + mu) * faces[Fpp].normale.x() * faces[Fpp].normale.z(),  (lambda + mu) * faces[Fpp].normale.x() * faces[Fpp].normale.y(),  (lambda + mu) * faces[Fpp].normale.y() * faces[Fpp].normale.y() + mu, (lambda + mu) * faces[Fpp].normale.y() * faces[Fpp].normale.z(), (lambda + mu) * faces[Fpp].normale.x() * faces[Fpp].normale.z(), (lambda + mu) * faces[Fpp].normale.y() * faces[Fpp].normale.z(), (lambda + mu) * faces[Fpp].normale.z() * faces[Fpp].normale.z() + mu;
    A_FppFpp *= faces[Fpp].S / V;

    //Assemblage de la matrice
    Mat.block<3,3>(0,0) = A_FF;
    Mat.block<3,3>(0,3) = A_FFp;
    Mat.block<3,3>(0,6) = A_FFpp;
    Mat.block<3,3>(3,0) = A_FpF;
    Mat.block<3,3>(3,3) = A_FpFp;
    Mat.block<3,3>(3,6) = A_FpFpp;
    Mat.block<3,3>(6,0) = A_FppF;
    Mat.block<3,3>(6,3) = A_FppFp;
    Mat.block<3,3>(6,6) = A_FppFpp;

    //Assemblage du second membre
    b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.),  ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,0.,1.),  ((-contrainte) * faces[Fpp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fpp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fpp].normale) * Vector_3(0.,0.,1.);

    //Modifs pour prendre en compte BC de Dirichlet s'il y en a
    if(faces[F].BC == 1) {
      Mat(2,0) = 0.;
      Mat(2,1) = 0.;
      Mat(2,2) = 1.;
      Mat(2,3) = 0.;
      Mat(2,4) = 0.;
      Mat(2,5) = 0.;
      Mat(2,6) = 0.;
      Mat(2,7) = 0.;
      Mat(2,8) = 0.;
      b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), displacement_BC_bis(faces[F].centre, solide[faces[F].voisins[0]].Dx, t, 0.),  ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,0.,1.), ((-contrainte) * faces[Fpp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fpp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fpp].normale) * Vector_3(0.,0.,1.);
    }
    else if(faces[Fp].BC == 1) {
      Mat(5,0) = 0.;
      Mat(5,1) = 0.;
      Mat(5,2) = 0.;
      Mat(5,3) = 0.;
      Mat(5,4) = 0.;
      Mat(5,5) = 1.;
      Mat(5,6) = 0.;
      Mat(5,7) = 0.;
      Mat(5,8) = 0.;
      b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.),  ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.),  displacement_BC_bis(faces[Fp].centre, solide[faces[Fp].voisins[0]].Dx, t, 0.), ((-contrainte) * faces[Fpp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fpp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fpp].normale) * Vector_3(0.,0.,1.);
    }
    else if(faces[Fpp].BC == 1) {
      Mat(8,0) = 0.;
      Mat(8,1) = 0.;
      Mat(8,2) = 0.;
      Mat(8,3) = 0.;
      Mat(8,4) = 0.;
      Mat(8,5) = 0.;
      Mat(8,6) = 0.;
      Mat(8,7) = 0.;
      Mat(8,8) = 1.;
      b << ((-contrainte) * faces[F].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[F].normale) * Vector_3(0.,0.,1.),  ((-contrainte) * faces[Fp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,1.,0.), ((-contrainte) * faces[Fp].normale) * Vector_3(0.,0.,1.), ((-contrainte) * faces[Fpp].normale) * Vector_3(1.,0.,0.), ((-contrainte) * faces[Fpp].normale) * Vector_3(0.,1.,0.), displacement_BC_bis(faces[Fp].centre, solide[faces[Fp].voisins[0]].Dx, t, 0.);
    }

    //Inversion du système !
    typedef Eigen::Matrix<double, 9, 9> Matrix9x9;
    Eigen::FullPivLU<Matrix9x9> lu(Mat);
    if( lu.rank() == 9) //Test voir si système inversible...
      x = Mat.lu().solve(b); 
    else { //Calcul de la pseudo-inverse pour minimisation de l'écart aux moindres carrés.
      Eigen::CompleteOrthogonalDecomposition<Matrix9x9> mat(Mat);
      x = mat.solve(b);
    }
    faces[F].I_Dx.vec[0] = x(0); faces[F].I_Dx.vec[1] = x(1); faces[F].I_Dx.vec[2] = x(2); //Première face de Neumann
    faces[Fp].I_Dx.vec[0] = x(3); faces[Fp].I_Dx.vec[1] = x(4); faces[Fp].I_Dx.vec[2] = x(5); //Deuxième face de Neumann
    faces[Fpp].I_Dx.vec[0] = x(6); faces[Fpp].I_Dx.vec[1] = x(7); faces[Fpp].I_Dx.vec[2] = x(8); //Deuxième face de Neumann
  }
}


void Solide::Forces_internes(const double& dt, const double& t, const double& T){ //Calcul des forces pour chaque particule
  stresses(t, T);
  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++) //Remet à zéro toutes les forces
    P->Fi = Vector_3(0.,0.,0.);
  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++){
    for(int i=0 ; i<P->faces.size() ; i++){
      int num_face = P->faces[i]; //Numéro de la face dans l'ensemble des faces contenu dans le solide
      /*int part_1 = faces[num_face].voisins[0];
	int part_2 = faces[num_face].voisins[1];*/
      if(faces[num_face].BC == 0) { //Forces sur faces internes
        int aux_1 = faces[num_face].reconstruction[0];
	double c_aux_1 = faces[num_face].c_reconstruction[0];
	int aux_2 = faces[num_face].reconstruction[1];
	double c_aux_2 = faces[num_face].c_reconstruction[1];
	int aux_3 = faces[num_face].reconstruction[2];
	double c_aux_3 = faces[num_face].c_reconstruction[2];
	int aux_4 = faces[num_face].reconstruction[3];
	double c_aux_4 = faces[num_face].c_reconstruction[3];
	//cout << "coords bary : " << c_aux_1 << " " << c_aux_2 << " " << c_aux_3 << " " << c_aux_4 << endl;
	
	//Sortir le sens de toutes les forces comme il faut...
	Vector_3 nIJ = faces[num_face].normale;
	/*int voisin;
	if(P->id == faces[num_face].voisins[0])
	  voisin = faces[num_face].voisins[1];
	else if(P->id == faces[num_face].voisins[1])
	voisin = faces[num_face].voisins[0];*/
	if(nIJ * Vector_3(P->x0, faces[num_face].centre) < 0.)
	  nIJ = -nIJ; //Normale pas dans le bon sens...

	//P->Fi = P->Fi +  2 * mu * (solide[voisin].Dx - P->Dx); //OK flux à 2 points
	//P->Fi = P->Fi + faces[num_face].S * (solide[part_1].contrainte + solide[part_2].contrainte) / 2. * nIJ; //OK Voronoi
	P->Fi = P->Fi + faces[num_face].S * P->contrainte * nIJ; //Force sur la particule sommée sur toutes les faces. Pour Tetra !
	solide[aux_1].Fi = solide[aux_1].Fi - faces[num_face].S * c_aux_1 * P->contrainte * nIJ; //Force sur autres particules en interaction avec P
	solide[aux_2].Fi = solide[aux_2].Fi - faces[num_face].S * c_aux_2 * P->contrainte * nIJ;
	solide[aux_3].Fi = solide[aux_3].Fi - faces[num_face].S * c_aux_3 * P->contrainte * nIJ;
	solide[aux_4].Fi = solide[aux_4].Fi - faces[num_face].S * c_aux_4 * P->contrainte * nIJ;
      }
      else if(faces[num_face].BC != 0) { //Calcul forces sur DDL bords Dirichlet. Vaut 0 exactement en Neumann Homogène. A vérifier... // == 1
	//cout << "Face au bord" << endl;
	int part = faces[num_face].voisins[0];
	Vector_3 nIJ = faces[num_face].normale;
	/*if((faces[num_face].S * solide[part].contrainte * nIJ).squared_length() > 1.)
	  cout << "Pb avec force sur bord de Neumann : " << (faces[num_face].S * solide[part].contrainte * nIJ).squared_length() << endl;*/
	P->Fi = P->Fi + faces[num_face].S * solide[part].contrainte * nIJ; //pow(10., 7.) * nIJ;
      }
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
      vtk << faces[*F].vertex.size();
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
    //vtk << 9 << endl; //Pour Quad
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
  //Tangent 1
  /*vtk << "VECTORS tan_1 double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
    vtk << faces[*F].vec_tangent_1 << endl;
  }
  vtk << "\n";
  //Tangent 2
  vtk << "VECTORS tan_2 double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
    vtk << faces[*F].vec_tangent_2 << endl;
  }
  vtk << "\n";*/
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
  /*vtk << "SCALARS p double 1" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      vtk << P->def_plas_cumulee << endl;
  }
  vtk << "\n";*/
  //Contrainte Von Mises
  vtk << "SCALARS Von_Mises double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
      vtk << (P->contrainte).VM() << endl;
    }
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
