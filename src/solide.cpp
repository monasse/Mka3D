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
//#include <eigen3/Eigen/LU> //Sert pour inversion systeÃ¨me linÃ©aire pour condition Neumann
#include <eigen3/Eigen/Dense>
#include <cmath>
#include <stdlib.h> //Pour utiliser system
#ifndef SOLIDE_CPP
#define SOLIDE_CPP

Solide::Solide(const double& E, const double& nu, const double& B1, const double& n1, const double& A1, const double& H1, const double& G, const int& recon){
  lambda = E * nu / (1.+nu) / (1. - 2.*nu);
  mu = E / 2. / (1.+nu);
  A = A1;
  B = B1;
  n = n1;
  H = H1;
  Gc = G;
  h = 0.;
  reconstruction = recon;
}

Solide::Solide(){
  lambda = 0.;
  mu = 0.;
  A = 0.;
  B = 0.;
  n = 0.;
  H = 0.;
  h = 0.;
  Gc = 0.;
  reconstruction = 0;
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

// void Solide::Init(const char* s1, const char* s2, const char* s3, const bool& rep, const int& numrep, const double& rho){ //pour Tetgen
//   std::ifstream noeuds(s1,ios::in);
//   std::ifstream elements(s2,ios::in);
//   //std::ifstream voisins(s3,ios::in);
//   std::ifstream import_faces(s3,ios::in);
//   if(not(noeuds) || not(elements) /*|| not(voisins)*/ || not(import_faces))
//     cout <<"ouverture du maillage ratee" << endl;

//   //Importation des vertex
//   string s;
//   char aux;
//   string ligne;
//   int nbr_vertex;
//   getline(noeuds, ligne); //Lecture de la premiÃ¨re ligne du fichier des vertex
//   istringstream  stdd(ligne);
//   stdd >> nbr_vertex; //Nombre de vertex

//   for(int i=0; i < nbr_vertex ; i++) {
//     getline(noeuds, ligne);
//     istringstream  stm(ligne);
//     int id; //NumÃ©ro du vertex
//     double x,y,z;
//     stm >> id >> x >> y >> z;
//     vertex.push_back(Vertex(Point_3(x,y,z), id)); //Vertex sont donnÃ©s dans l'ordre
//   }

//   //Importation des Particules
//   getline(elements, ligne);
//   istringstream stde(ligne);
//   int nbr_elements;
//   stde >> nbr_elements; //On stocke le nombre d'Ã©lÃ©ments
//   for(int i=0; i < nbr_elements ; i++) {
//     getline(elements, ligne);
//     istringstream  stm(ligne);
//     int id;
//     stm >> id;
//     int v1,v2,v3,v4;
//     stm >> v1 >> v2 >> v3 >> v4;
//     //Ajout des vertex de la particule
//     Particule p(id);
//     p.vertices.push_back(v1);
//     p.vertices.push_back(v2);
//     p.vertices.push_back(v3);
//     p.vertices.push_back(v4);

//     //Calcul des quantitÃ©s volumiques (liÃ©es particule)
//     p.barycentre(this, 4); //Calcul du barycentre
//     p.volume(this, 4); //calcul du volume
//     p.m = rho * p.V;

//     //Ajout de la particule dans le solide
//     solide.push_back(p);
//   }
  
//   //Importation des faces et des connectivitÃ©s
//   getline(import_faces, ligne);
//   istringstream stdf(ligne);
//   int nbr_faces;
//   stdf >> nbr_faces;
//   for(int i=0; i<nbr_faces ; i++) {
//     getline(import_faces, ligne);
//     istringstream  stm(ligne);
//     int id;
//     int v1,v2,v3; //Les vertex de la face
//     int BC; //Condition de bord
//     int part_1, part_2; //Le numÃ©ro des particules voisines
//     stm >> id >> v1 >> v2 >> v3 >> BC >> part_1 >> part_2; //NumÃ©ro de la face, des vertex + 1 et des voisins (bon numÃ©ro)
//     Face f;
//     f.id = id;
//     f.vertex.push_back(v1); //Ajout du numÃ©ro des vertex
//     f.vertex.push_back(v2);
//     f.vertex.push_back(v3);
//     f.BC = BC;
//     f.type = 2; //Triangle pour tegen
//     /*if(part_1 == -1 || part_2 == -1)
//       f.BC = -1; //face au bord
//     else {
//       f.BC = 0; //Face pas au bord
//       f.D0 = sqrt(Vector_3(solide[part_1].x0, solide[part_2].x0).squared_length());
//       }*/
//     if(part_1 >=0) { //Ajout du numÃ©ro des voisins dans la face
//       f.voisins.push_back(part_1);
//       f.voisins.push_back(part_2);
//     }
//     else { //La premiÃ¨re valeur de voisin est la seule particule qui contient la face sur le bord
//       f.voisins.push_back(part_2);
//       f.voisins.push_back(part_1);
//     }
//     //bool calcul_normales = false;
//     if(part_1 >= 0) {
//       solide[part_1].faces.push_back(f.id); //Ajout du numÃ©ro de la face dans la liste ds voisins de chaque particule
//       if(f.BC != 0) {
// 	if(f.BC == 1)
// 	  solide[part_1].BC = f.BC;
// 	else if(solide[part_1].BC <= 0)
// 	  solide[part_1].BC = f.BC;
//       }
//     }
//     if(part_2 >= 0) {
//       solide[part_2].faces.push_back(f.id); //Ajout du numÃ©ro de la face dans la liste ds voisins de chaque particule
//       if(f.BC != 0) {
// 	if(f.BC == 1)
// 	  solide[part_1].BC = f.BC;
// 	else if(solide[part_1].BC <= 0)
// 	  solide[part_1].BC = f.BC;
//       }
//     }    
//     f.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face

//     //VÃ©rification du sens de la normale
//     if(part_1 != -1 && part_2 != -1) { //Face pas au bord
//       if(Vector_3(solide[part_1].x0, solide[part_2].x0) * f.normale  < 0.)
// 	f.normale = -f.normale;
//     }
//     if(part_1 == -1 || part_2 == -1) { //Face au bord
//       Vector_3 bonne_direction = Vector_3(solide[f.voisins[0]].x0, f.centre);
//       if(bonne_direction * f.normale < 0.)
// 	f.normale = -f.normale;
//     }
//     faces.push_back(f);
//   }

//   //Calcul du tetrahÃ¨dre associÃ© Ã  chaque face pour le calcul du gradient
//   for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){ //Boucle sur toutes les faces
//     if(F->BC == 0 && not(F->split)) {
//       //cout << "Face : " << F->id << endl;
//       bool test = voisins_face(F->id);
//       if(not(test)) {
// 	cout << "Face : " << F->id << " Pas de tetra associe a une face" << endl;
// 	throw std::invalid_argument( "Pas de tetra associe a une face" );
// 	/*cout << "Centre Face : " << F->centre << endl;
// 	cout << "Barycentre Voisin A : " << solide[F->voisins[0]].x0 << endl;
// 	cout << "Barycentre Voisin A : " << solide[F->voisins[0]].x0 << endl;
// 	cout << "Barycentre Voisin B : " << solide[F->voisins[1]].x0 << endl;*/
//       }
//     }
//   }
// }

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
    vertex.push_back(Vertex(Point_3(x,y,z), id-1)); //Vertex sont donnÃ©s dans l'ordre
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
      //cylindre
      //if(tag_2 == 11 || tag_2 == 33) //Dirichlet
      if(tag_1 == 40) //Traction 1
	F.BC = 1;
      else if(tag_1 == 43) //Neumann homogène
	F.BC = -1;
      else if(tag_1 == 39) { //Fissure déjà présente et Neumann Homogène
	F.BC = -2;
	F.fissure = true;
      }
      else if(tag_1 == 41) //Traction 2
	F.BC = 2;
      else if(tag_1 == 42) //Déplacements dans le plan bloqués
	F.BC = 3;
      //Poutre section carré
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
      vertex[v1-1].particules.push_back(solide.size());
      vertex[v2-1].particules.push_back(solide.size());
      vertex[v3-1].particules.push_back(solide.size());
      vertex[v4-1].particules.push_back(solide.size());

      //Calcul des quantités volumiques (liées particules)
      p.barycentre(this, type); //Calcul du barycentre
      p.volume(this, type); //calcul du volume
      p.m = rho * p.V;
      p.calcul_diametre(this); //Calcul du diamètre de la particule

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

      //Calcul des quantitÃ©s volumiques (liÃ©es particule)
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

  taille_maillage(); //Calcul de la taille du maillage

  //cout << "nombre particules : " << solide.size() << endl;
  //cout << "Nombre total de faces : " << faces.size() << endl;

  //Création des connectivités entre éléments
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){ //Boucle sur toutes les faces
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      if(P->contient_face(*F)) { //On ajoute les numÃ©ros de la face dans la particule et rÃ©ciproquement
	//cout << "Particule : " << P->id << endl;
	P->faces.push_back(F->id);
	(F->voisins).push_back(P->id);
	if(F->BC != 0 && F->BC != -2)
	  (F->voisins).push_back(-1); //Car face au bord
	if(type == 5)
	  (F->c_reconstruction).push_back(0.5); //Maillage de Voronoi
	if(F->voisins.size() > 2) {
	  cout << "Numero face : " << F->id << endl;
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

  //On fait en sorte que la normale aille de voisins[0] vers voisins[1]
  /*for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
    if( F->normale * Vector_3(solide[F->voisins[0]].x0, solide[F->voisins[1]].x0) < 0.)
      F->normale = -F->normale;
      }*/

  cout << "Debut masses faces" << endl;
  //Boucle pour donner une masse aux faces
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
    if(F->BC != 0) {
      int voisin = F->voisins[0];
      F->m = sqrt(pow(F->S * (F->centre - solide[voisin].x0) * F->normale / 3. * rho, 2.));
      //cout << "Face : " << F->id << " masse : " << F->m << endl;
    }
    else {
      int voisin1 = F->voisins[0];
      int voisin2 = F->voisins[1];
      F->m = sqrt(pow(F->S * (F->centre - solide[voisin1].x0) * F->normale / 3. * rho + F->S * (F->centre - solide[voisin2].x0) * F->normale / 3. * rho, 2.));
      //cout << "Face : " << F->id << " masse : " << F->m << endl;
      F->masses.push_back( sqrt(pow(F->S * (F->centre - solide[voisin1].x0) * F->normale / 3. * rho, 2.)) );
      F->masses.push_back( sqrt(pow(F->S * (F->centre - solide[voisin2].x0) * F->normale / 3. * rho, 2.)) );
    }
    cout << "ok face : " << F->id << endl; 
  }

  //Il faudra refaire un peu les connectivités mais ça va...
  //Faire ici une boucle sur les particules pour vérifier le nombre de faces de Neumann qu'elles ont. S'il y en a plus d'une splitter la particule et avoir un hanging node. Les tetras voisins auront des hanging nodes mais c'est pas trop grave
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    int nb_faces = 0; //Sert à récupérer le numéro des faces avec BC de Neumann si nécessaire
    for(int i=0 ; i < P->faces.size() ; i++){
      int f = P->faces[i];
      //cout << "Num part : " << P->id << " Num faces : " << f << endl;
      if(faces[f].BC != 0) //car on ne fait ces calculs sur toutes les faces au bord
	nb_faces++;
    }
    if(nb_faces > 1) { //On appelle la fonction de splitting
      //cout << "On appelle le splitting !" << endl;
      //splitting_elements(P->id, rho);
    }
  }

  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      for(int i=0 ; i < P->faces.size() ; i++){
	int f = P->faces[i];
	if(faces[f].split){
	  cout << "Gros probleme !" << endl;
	  getchar();
	}
      }
    }
  }
  //cout << "On a splitte tout ce qu'il faut !" << endl;
  //getchar();

  //Verification de la connectivite
  //Verification de la numérotation des particules
  for(int i=0;i<solide.size();i++){
    const Particule& P = solide[i];
    if(i!=P.id){
      cout << "probleme numerotation particule i=" << i << " P.id=" << P.id << endl;
      getchar();
    }
  }
  //Verification de la numérotation des faces
  for(int i=0;i<faces.size();i++){
    const Face& F = faces[i];
    if(i!=F.id){
      cout << "probleme numerotation face i=" << i << " F.id=" << F.id << endl;
      getchar();
    }
  }
  //Verification de la numérotation des vertex
  for(int i=0;i<vertex.size();i++){
    const Vertex& V = vertex[i];
    if(i!=V.num){
      cout << "probleme numerotation vertex i=" << i << " V.num=" << V.num << endl;
      getchar();
    }
  }
  //Verification de la connectivite faces/particules
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    //if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	if(faces[*F].voisins[0]!=P->id && faces[*F].voisins[1]!=P->id){
	  cout << "Probleme connectivite face=" << *F << " voisins=" << faces[*F].voisins[0] << " " << faces[*F].voisins[1] << " part.id=" << P->id << endl;
	  getchar();
	}
      }
      //}
  }
  //Verification de la connectivite particules/faces
  for(int i=0;i<faces.size();i++){
    Face F = faces[i];
    if(F.voisins.size()!=2){
      cout << "La face " << i << " n'a pas le bon nombre de voisins ! size=" << F.voisins.size() << endl;
      for(std::vector<int>::iterator P=F.voisins.begin();P!=F.voisins.end();P++){
	cout << *P << " ";
      }
      cout << endl;
      getchar();
    }
  }
  for(int i=0;i<faces.size();i++){
    Face F = faces[i];
    if(F.voisins[0]>=0){
      Particule P = solide[F.voisins[0]];
      bool test=true;
      for(std::vector<int>::iterator f=P.faces.begin();f!=P.faces.end() && test;f++){
	test = (*f!=i);
      }
      if(test){
	cout << "Probleme connectivite face=" << i << " voisins[0]=" << F.voisins[0] << " P.id=" << P.id << endl;
      }
    }
    if(F.voisins[1]>=0){
      Particule P = solide[F.voisins[1]];
      bool test=true;
      for(std::vector<int>::iterator f=P.faces.begin();f!=P.faces.end() && test;f++){
	test = (*f!=i);
      }
      if(test){
	cout << "Probleme connectivite face=" << i << " voisins[1]=" << F.voisins[1] << " P.id=" << P.id << endl;
      }
    }
  }
  //Verification de la connectivite faces/vertex
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
      for(std::vector<int>::iterator V=faces[*F].vertex.begin();V!=faces[*F].vertex.end();V++){
	bool test=true;
	for(std::vector<int>::iterator v=P->vertices.begin();v!=P->vertices.end() && test;v++){
	  test = (*v!=*V);
	}
	if(test){
	  cout << "Probleme connectivite face=" << *F << " vertex=" << *V << " particule=" << P->id << endl;
	  cout << "Points de la face ";
	  for(std::vector<int>::iterator Vt=faces[*F].vertex.begin();Vt!=faces[*F].vertex.end();Vt++){
	    cout << *Vt << " ";
	  }
	  cout << endl;
	  cout << "Points de la particule ";
	  for(std::vector<int>::iterator v=P->vertices.begin();v!=P->vertices.end() && test;v++){
	    cout << *v << " ";
	  }
	  cout << endl;
	  getchar();
	}
      }
    }
    for(std::vector<int>::iterator v=P->vertices.begin();v!=P->vertices.end();v++){
      bool test=true;
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end() && test;F++){
	for(std::vector<int>::iterator V=faces[*F].vertex.begin();V!=faces[*F].vertex.end() && test;V++){
	  test = (*v!=*V);
	}
      }
      if(test){
	cout << "Probleme connectivite vertex=" << *v << " particule=" << P->id << endl;
	getchar();
      }
    }
  }
  //Verification de la connectivite particules/vertex
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator V=P->vertices.begin();V!=P->vertices.end();V++){
      Vertex v=vertex[*V];
      bool test=true;
      for(std::vector<int>::iterator p=v.particules.begin();p!=v.particules.end() && test;p++){
	test = (*p!=P->id);
      }
      if(test){
	cout << "Probleme connectivite particule=" << P->id << " vertex=" << *V << endl;
	for(std::vector<int>::iterator p=v.particules.begin();p!=v.particules.end() && test;p++){
	  cout << *p << " " ;
	}
	cout << endl;
	getchar();
      }
    }
  }
  //Verification de la connectivite vertex/particules
  for(std::vector<Vertex>::iterator v=vertex.begin();v!=vertex.end();v++){
    for(std::vector<int>::iterator p=v->particules.begin();p!=v->particules.end();p++){
      Particule P = solide[*p];
      bool test=true;
      for(std::vector<int>::iterator V=P.vertices.begin();V!=P.vertices.end() && test;V++){
	test = (*V!=v->num);
      }
      if(test){
	cout << "Probleme connectivite vertex=" << v->num << " particule=" << *p << " P.id=" << P.id << endl;
      }
    }
  }
  
  cout << "Connectivite OK" << endl;
  
  
  
  
  
  
  //Calcul du tetrahèdre associé à chaque face pour le calcul du gradient
  /*if(type == 4) { //Seulement pour tetra
    for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){ //Boucle sur toutes les faces
      if(F->BC == 0 && not(F->split)) {
	//cout << "Face : " << F->id << endl;
	bool test = voisins_face(F->id);
	if(not(test)) {
	  cout << "Face : " << F->id << " Extrapolation sur la face avec voisins_face" << endl;
	  //throw std::invalid_argument( "Pas de tetra associe a une face" );
	}
      }
    }
    }*/

  /*if(type == 4) { //Utilisation du Delaunay pour trouver les tÃ©tras associÃ©s Ã  chaque face
    bool calcul_pre_traitement = false;
    std::ifstream pre_traitement("tetra_faces.txt"); //,ios::in);

    bool test = getline(pre_traitement, ligne);
    //cout << "test : " << test << endl;
    //Si pas de pretraitement
    //if(test) {
    if(false){ //Obligation de prétraitement !
      bool aux = true;
      for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
	if(aux && F->BC == 0 && not(F->split)) { //Importation des tetras associés aux faces internes
	  istringstream  stm(ligne);
	  int id, ele1,ele2,ele3,ele4;
	  double c1,c2,c3,c4;
	  stm >> id >> ele1 >> ele2 >> ele3 >> ele4 >> c1 >> c2 >> c3 >> c4;
	  if(id != F->id){
	    cout << "probleme numerotation face id=" << id << " F.id=" << F->id << endl;
	    getchar();
	  }
	  F->reconstruction.clear();
	  F->reconstruction.push_back(ele1);
	  F->reconstruction.push_back(ele2);
	  F->reconstruction.push_back(ele3);
	  F->reconstruction.push_back(ele4);
	  if(solide[ele1].split){
	    cout << "ele1 split !" << endl;
	    getchar();
	  }
	  if(solide[ele2].split){
	    cout << "ele2 split !" << endl;
	    getchar();
	  }
	  if(solide[ele3].split){
	    cout << "ele3 split !" << endl;
	    getchar();
	  }
	  if(solide[ele4].split){
	    cout << "ele4 split !" << endl;
	    getchar();
	  }
	  F->c_reconstruction.clear();
	  F->c_reconstruction.push_back(c1);
	  F->c_reconstruction.push_back(c2);
	  F->c_reconstruction.push_back(c3);
	  F->c_reconstruction.push_back(c4);
	  //cout << "face : " << F->id << endl;

	  aux = getline(pre_traitement, ligne);
	}
      }
    }
    else {
      //On génère le fichier pour qhull qui va permettre ed calculer le delaunay avec barycentre des mailles
      std::ofstream centres("cell_centres.txt",ios::out);
      centres << 3 << " Cell-centres" << endl; //Pour indiquer dimension des points
      int nb_centres = 0;
      for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
	if(not(P->split)) //Particule pas splitée
	  nb_centres++; //On va mettre la particule dans le fichier
      }
      centres << nb_centres << endl;
      for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
	if(not(P->split)) //Particule pas splitée
	  centres << P->x0 << endl; //On met la particule dans le fichier
      }

      //On va appeler la fonction qhull voulu sur le fichier précédent avec commande console
      int rez = system("qdelaunay i Qt < cell_centres.txt > delaunay.txt"); //Commande bash pour charger les cell-centres

      //On importe la tetrahedrisation de Delaunay
      std::ifstream delaunay("delaunay.txt",ios::in);
      if(not(delaunay))
	throw std::invalid_argument( "Pas de tetrahedrisation de Delaunay !" );
      getline(delaunay, ligne); //On enlève la ligne avec le nombre d'éléments
      std::vector<Particule> tetra_delau; //Ensemble des particules contenant la tetra
      while(getline(delaunay, ligne)) { //Importation de la tetraedrisation de Delaunay
	istringstream  stm(ligne);
	int ele1,ele2,ele3,ele4;
	stm >> ele1 >> ele2 >> ele3 >> ele4;
	Particule P;
	P.vertices.push_back(ele1); //NumÃ©ros des Elements du solide qui forment chacun des tetra
	P.vertices.push_back(ele2);
	P.vertices.push_back(ele3); // - 1
	P.vertices.push_back(ele4);
	tetra_delau.push_back(P);
      }
      //cout << "ok Delaunay" << endl;

      //On sort dans un fichier les éléments associés à chaque face ainsi que les coodonnées barycentriques correspondantes
      std::ofstream tetra_faces("tetra_faces.txt",ios::out); //Sortie pour stocker tetras associés aux faces
      std::ofstream face_pb("face_pb.txt",ios::out); //Sorties pour les faces qui pose pb
      //Recherche du tetraÃ¨dre associÃ© Ã  chaque face
      for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
	bool tetra_ok = false;
	if(F->BC == 0 && not(F->split)) {
	  //cout << "Voisins : " << F->voisins[0] << " " << F->voisins[1] << endl;
	  // for(std::vector<Particule>::iterator P=tetra_delau.begin();P!=tetra_delau.end();P++){
	  //   int part_1 = P->vertices[0];
	  //   int part_2 = P->vertices[1];
	  //   int voisin1 = P->vertices[2];
	  //   int voisin2 = P->vertices[3];
	  
	  //   double c1 = (Vector_3(solide[part_2].x0, F->centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_2].x0, solide[part_1].x0) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0) ));
	  //   double c2 = (Vector_3(solide[part_1].x0, F->centre) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_1].x0, solide[part_2].x0) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0) ));
	  //   double c3 = (Vector_3(solide[part_2].x0, F->centre) * cross_product(Vector_3(solide[part_2].x0, solide[part_1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0)) ) / (Vector_3(solide[part_2].x0, solide[voisin1].x0) * cross_product(Vector_3(solide[part_2].x0, solide[part_1].x0), Vector_3(solide[part_2].x0, solide[voisin2].x0) ));
	  //   double c4 = (Vector_3(solide[part_2].x0, F->centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[part_1].x0)) ) / (Vector_3(solide[part_2].x0, solide[voisin2].x0) * cross_product(Vector_3(solide[part_2].x0, solide[voisin1].x0), Vector_3(solide[part_2].x0, solide[part_1].x0) ));

	  //   if( c1 >= 0. && c2 >= 0. && c3 >= 0. && c4 >= 0. && c1 < 1. && c2 < 1. && c3 < 1. && c4 < 1.) {
	  //     F->reconstruction.push_back(part_1);
	  //     F->reconstruction.push_back(part_2);
	  //     F->reconstruction.push_back(voisin1);
	  //     F->reconstruction.push_back(voisin2);
	  //     F->c_reconstruction.push_back(c1);
	  //     F->c_reconstruction.push_back(c2);
	  //     F->c_reconstruction.push_back(c3);
	  //     F->c_reconstruction.push_back(c4);
	  //     //cout << F->id << endl;
	  //     //if(c1 >= 0.&& c2 >= 0. && c3 >= 0.)
	  //     //cout << c1 << " " << c2 << " " << c3 << " " << c4 << " " << c1 + c2 + c3 + c4 - 1. << endl;
	  //     //cout << "face : " << F->id << " ok !" << endl;
	  //     tetra_ok = true;
	  //     //cout << "Face : " << F->id << " ok !" << endl;
	  //     break;
	  //   }
	  // }
	  // if( not(tetra_ok)) {
	  if(true){
	    //cout << "Face : " << F->id << endl;
	    //throw std::invalid_argument( " pas de tetra associe a une face !" );
	    face_pb << F->id << " " << vertex[F->vertex[0]].pos << " " << vertex[F->vertex[1]].pos << " " << vertex[F->vertex[2]].pos << endl;
	    //throw std::invalid_argument( " pas de tetra associe a la face !" );
	    bool test = voisins_face(F->id); //Dans ce cas, on va faire de l'extrapolation et utiliser l'ancienne méthode...
	    //cout << "face ok !" << endl;
	    if(not(test)) {
	      F->face_pb = true;
	      cout << "Nbr Faces : " << faces.size() << " Face : " << F->id << endl;
	      //getchar();
	      //throw std::invalid_argument( "Pas de tetra associe a une face" );
	      cout << "Pas de tetra associe a une face" << endl;
	    }
	  }
      
	//Sortie contenant les Ã©lÃ©ments associÃ©s Ã  chaque face pour faire le prÃ©-traitement une seule fois
	//cout << "Enregistrement pret !" << endl;
	  tetra_faces << F->id << " " << F->reconstruction[0] << " " << F->reconstruction[1] << " " << F->reconstruction[2] << " " << F->reconstruction[3] << " " << F->c_reconstruction[0] << " " << F->c_reconstruction[1] << " " << F->c_reconstruction[2] << " " << F->c_reconstruction[3] << endl;
	}
      }
    }
  }*/
}

void Solide::splitting_elements(const int& num_part, const double& rho) {
  Particule* P = &(solide[num_part]);
  std::vector<int> out; //Face neumann
  std::vector<int> in; //Inner faces

  for(int i=0 ; i < P->faces.size() ; i++){
    int f = P->faces[i];
    if(faces[f].BC == 0)
      in.push_back(f);
    else
      out.push_back(f);
  }
  //cout << "Taille in : " << in.size() << " taille out : " << out.size() << endl;


  //Trouver edge commun aux deux faces internes
  std::vector<int> common_edge;
  for(int i=0; i < faces[in[0]].vertex.size() ; i++) {
    for(int j=0; j < faces[in[1]].vertex.size() ; j++) {
      if(faces[in[0]].vertex[i] == faces[in[1]].vertex[j])
	common_edge.push_back(faces[in[1]].vertex[j]);
    }
  }
  if(common_edge.size() > 2 || common_edge.size() < 2)
    throw std::invalid_argument( "Pas d'edge commun entre les 2 faces a splitter" );
  //Ordonnancement de common_edge selon l'ordre de out
  bool test=true;
  for(std::vector<int>::iterator v=faces[out[0]].vertex.begin();v!=faces[out[0]].vertex.end() && test;v++){
    test = (*v!=common_edge[0]);
  }
  if(test){
    int tmp = common_edge[0];
    common_edge[0] = common_edge[1];
    common_edge[1] = tmp;
  }
  for(std::vector<int>::iterator v=faces[out[0]].vertex.begin();v!=faces[out[0]].vertex.end() && test;v++){
    test = (*v!=common_edge[0]);
  }
  if(test){
    cout << "la correction n'a pas fontionne common_edge[0]=" << common_edge[0] << endl;
    getchar();
  }
  for(std::vector<int>::iterator v=faces[out[1]].vertex.begin();v!=faces[out[1]].vertex.end() && test;v++){
    test = (*v!=common_edge[1]);
  }
  if(test){
    cout << "la correction n'a pas fontionne common_edge[1]=" << common_edge[1] << endl;
    getchar();
  }
  

  //On calcule la moitiéde l'edge et on split les faces
  double id = vertex.size(); //Numéro du vertex qu'on va ajouter
  Vertex demi_edge = Vertex(Point_3(0.5 * (vertex[common_edge[0]].pos.x() + vertex[common_edge[1]].pos.x()),0.5 * (vertex[common_edge[0]].pos.y() + vertex[common_edge[1]].pos.y()), 0.5 * (vertex[common_edge[0]].pos.z() + vertex[common_edge[1]].pos.z()) ), id); //Nouveau vertex
  vertex.push_back(demi_edge); //Vertex sont donnés dans l'ordre

  //on stocke les vertex qui ne sont pas sur l'edge splité
  std::vector<int> vertex_common_part_out, vertex_part_2;
  // for(int i=0; i < faces[out[0]].vertex.size() ; i++) {
  //   if(faces[out[0]].vertex[i] != common_edge[0])
  //     vertex_common_part_out.push_back(faces[out[0]].vertex[i]);
  // }
  for(int i=0; i < faces[out[0]].vertex.size() ; i++) {
    for(int j=0; j < faces[out[1]].vertex.size() ; j++) {
      if(faces[out[0]].vertex[i] == faces[out[1]].vertex[j])
	vertex_common_part_out.push_back(faces[out[1]].vertex[j]);
    }
  }

  //Création des 2 nouvelles particules issues du splitting
  Particule part_1, part_2;
  part_1.id = solide.size();
  part_2.id = solide.size() + 1;
  //Pas besoin de mettre les faces dans les particules. C'est fait plus tard dans l'algo. Il suffitd'y mettre les vertex
  part_1.vertices.push_back(id); //Ajout du nouvelle edge dans les 2 particules
  part_1.vertices.push_back(faces[out[0]].vertex[0]);
  part_1.vertices.push_back(faces[out[0]].vertex[1]);
  part_1.vertices.push_back(faces[out[0]].vertex[2]);
  //cout << "num_part splitte=" << num_part << endl;
  //cout << "part_1: id=" << id << " out[0]=" << faces[out[0]].vertex[0] << " " << faces[out[0]].vertex[1] << " " << faces[out[0]].vertex[2] << endl;
  part_2.vertices.push_back(id); //Ajout du nouvelle edge dans les 2 particules
  part_2.vertices.push_back(faces[out[1]].vertex[0]);
  part_2.vertices.push_back(faces[out[1]].vertex[1]);
  part_2.vertices.push_back(faces[out[1]].vertex[2]);
  //cout << "part_2: id=" << id << " out[1]=" << faces[out[1]].vertex[0] << " " << faces[out[1]].vertex[1] << " " << faces[out[1]].vertex[2] << endl;

  //Nouvelle face
  Face new_face;
  new_face.vertex.push_back(id); //Numéro du vertex issu du demi-edge
  new_face.vertex.push_back(vertex_common_part_out[0]);
  new_face.vertex.push_back(vertex_common_part_out[1]);
  new_face.BC = 0; //Face dans bulk du coup
  new_face.type = 2;
  new_face.id = faces.size(); ////Ajout de la face dans l'ensemble des faces du Solide
  new_face.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
  //cout << "New face ok" << endl;
  new_face.voisins.push_back(part_1.id);
  new_face.voisins.push_back(part_2.id);
  faces.push_back(new_face);
  
  //Faces Splitées (4 après splitting)
  //in[0]
  Face face1;
  if(faces[in[0]].vertex[0] != common_edge[0] && faces[in[0]].vertex[0] != common_edge[1])
    face1.vertex.push_back(faces[in[0]].vertex[0]); //Il fait choisir les vertex pas alignés avec le nouvel edge créé
  else if(faces[in[0]].vertex[1] != common_edge[0] && faces[in[0]].vertex[1] != common_edge[1])
    face1.vertex.push_back(faces[in[0]].vertex[1]);
  else if(faces[in[0]].vertex[2] != common_edge[0] && faces[in[0]].vertex[2] != common_edge[1])
    face1.vertex.push_back(faces[in[0]].vertex[2]);
  face1.vertex.push_back(common_edge[0]);
  face1.vertex.push_back(id);
  face1.type = 2;
  face1.BC = 0;
  face1.id = faces.size();
  face1.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
  //cout << "Face 1 ok" << endl;
  if(faces[in[0]].voisins[0] == num_part)
    face1.voisins.push_back(faces[in[0]].voisins[1]);
  else if(faces[in[0]].voisins[1] == num_part)
    face1.voisins.push_back(faces[in[0]].voisins[0]);
  face1.voisins.push_back(part_1.id);
  //cout << "Voisins : " << face1.voisins[0] << " " << face1.voisins[1] << endl;
  faces.push_back(face1);

  Face face2;
  if(faces[in[0]].vertex[0] != common_edge[0] && faces[in[0]].vertex[0] != common_edge[1])
    face2.vertex.push_back(faces[in[0]].vertex[0]);
  else if(faces[in[0]].vertex[1] != common_edge[0] && faces[in[0]].vertex[1] != common_edge[1])
    face2.vertex.push_back(faces[in[0]].vertex[1]);
  else if(faces[in[0]].vertex[2] != common_edge[0] && faces[in[0]].vertex[2] != common_edge[1])
    face2.vertex.push_back(faces[in[0]].vertex[2]);
  face2.vertex.push_back(common_edge[1]);
  face2.vertex.push_back(id);
  face2.type = 2;
  face2.BC = 0;
  face2.id = faces.size();
  face2.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
  //cout << "Face 2 ok" << endl;
  face2.voisins.push_back(face1.voisins[0]);
  face2.voisins.push_back(part_2.id);
  //cout << "Voisins : " << face2.voisins[0] << " " << face2.voisins[1] << endl;
  faces.push_back(face2);

  //in[1]
  Face face3;
  if(faces[in[1]].vertex[0] != common_edge[0] && faces[in[1]].vertex[0] != common_edge[1])
    face3.vertex.push_back(faces[in[1]].vertex[0]); //Il fait choisir les vertex pas alignés avec le nouvel edge créé
  else if(faces[in[1]].vertex[1] != common_edge[0] && faces[in[1]].vertex[1] != common_edge[1])
    face3.vertex.push_back(faces[in[1]].vertex[1]);
  else if(faces[in[1]].vertex[2] != common_edge[0] && faces[in[1]].vertex[2] != common_edge[1])
    face3.vertex.push_back(faces[in[1]].vertex[2]);
  face3.vertex.push_back(common_edge[0]);
  face3.vertex.push_back(id);
  face3.type = 2;
  face3.BC = 0;
  face3.id = faces.size();
  face3.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
  //cout << "Face 3 ok" << endl;
  if(faces[in[1]].voisins[0] == num_part)
    face3.voisins.push_back(faces[in[1]].voisins[1]);
  else if(faces[in[1]].voisins[1] == num_part)
    face3.voisins.push_back(faces[in[1]].voisins[0]);
  face3.voisins.push_back(part_1.id);
  //cout << "Voisins : " << face3.voisins[0] << " " << face3.voisins[1] << endl;
  faces.push_back(face3);

  Face face4;
  if(faces[in[1]].vertex[0] != common_edge[0] && faces[in[1]].vertex[0] != common_edge[1])
    face4.vertex.push_back(faces[in[1]].vertex[0]);
  else if(faces[in[1]].vertex[1] != common_edge[0] && faces[in[1]].vertex[1] != common_edge[1])
    face4.vertex.push_back(faces[in[1]].vertex[1]);
  else if(faces[in[1]].vertex[2] != common_edge[0] && faces[in[1]].vertex[2] != common_edge[1])
    face4.vertex.push_back(faces[in[1]].vertex[2]);
  face4.vertex.push_back(common_edge[1]);
  face4.vertex.push_back(id);
  face4.type = 2;
  face4.BC = 0;
  face4.id = faces.size();
  face4.comp_quantities(this); //Calcul de la normale sortante, surface et barycentre face
  //cout << "Face 4 ok" << endl;
  face4.voisins.push_back(face3.voisins[0]);
  face4.voisins.push_back(part_2.id);
  //cout << "Voisins : " << face4.voisins[0] << " " << face4.voisins[1] << endl;
  faces.push_back(face4);

  //Ajout des faces dans les 2 particules
  faces[out[0]].voisins[0] = part_1.id;
  part_1.faces.push_back(out[0]);
  part_1.faces.push_back(new_face.id);
  part_1.faces.push_back(face1.id);
  part_1.faces.push_back(face3.id);
  part_1.barycentre(this, 4); //Calcul du barycentre
  part_1.volume(this, 4); //calcul du volume
  part_1.m = rho * part_1.V;
  part_1.from_splitting = true;
  solide.push_back(part_1);

  faces[out[1]].voisins[0] = part_2.id;
  part_2.faces.push_back(out[1]);
  part_2.faces.push_back(new_face.id);
  part_2.faces.push_back(face2.id);
  part_2.faces.push_back(face4.id);
  part_2.barycentre(this, 4); //Calcul du barycentre
  part_2.volume(this, 4); //calcul du volume
  part_2.m = rho * part_2.V;
  part_2.from_splitting = true;
  solide.push_back(part_2);

  //Ajout des nouvelles faces dans particules environnantes pas splitées et retrait des faces splitées des particules
  int part_voisine_1 = face1.voisins[0];
  int part_voisine_2 = face3.voisins[0];
  //Première particule
  for(int i=0; i < solide[part_voisine_1].faces.size() ; i++) {
    int F = solide[part_voisine_1].faces[i];
    if(F == in[0]) {
      solide[part_voisine_1].faces[i] = face1.id; //On retire la face splitée et on met à la place une des 2 nouvelles
      break;
    }
  }
  solide[part_voisine_1].faces.push_back(face2.id);
  solide[part_voisine_1].vertices.push_back(id); //On ajoute le vertex créé au milieu de l'edge
  solide[part_voisine_1].impact_splitting = true;
  

  //Deuxième particule
  for(int i=0; i < solide[part_voisine_2].faces.size() ; i++) {
    int F = solide[part_voisine_2].faces[i];
    if(F == in[1]) {
      solide[part_voisine_2].faces[i] = face3.id; //On retire la face splitée et on met à la place une des 2 nouvelles
      break;
    }
  }
  solide[part_voisine_2].faces.push_back(face4.id);
  solide[part_voisine_2].vertices.push_back(id); //On ajoute le vertex créé au milieu de l'edge
  solide[part_voisine_2].impact_splitting = true;
  
  //Il reste à détruire (ou vider ?) la particule splittée ainsi que les 2 faces splittées !
  //Comment faire ? Mettre marqueur pour particule jetée à pas prendre en compte dans calcul des reconstructions etc ???
  //On indique que les faces et particules splitées ne doivent plus être prises en compte dans les calculs
  faces[in[0]].split = true;
  faces[in[1]].split = true;
  faces[in[0]].voisins[0] = -1;
  faces[in[0]].voisins[1] = -1;
  faces[in[0]].vertex.clear();
  faces[in[0]].reconstruction.clear();
  faces[in[0]].c_reconstruction.clear();
  faces[in[1]].voisins[0] = -1;
  faces[in[1]].voisins[1] = -1;
  faces[in[1]].vertex.clear();
  faces[in[1]].reconstruction.clear();
  faces[in[1]].c_reconstruction.clear();
  solide[num_part].split = true;
  solide[num_part].faces.clear();
  solide[num_part].vertices.clear();
  
  
  //getchar();
  //cout << "Particules crees : " << part_1.id << " " << part_2.id << endl;
  //cout << "Faces crees : " << new_face.id << " " << face1.id << " " << face2.id << " " << face3.id << " " << face4.id << endl;
  
}

bool Solide::face_existe(Face f) { //Renvoie vraie si la face testée est déjà dans faces
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
    if(*F == f) {
      return true;
    }
  }
  return false;
}

void Solide::taille_maillage() {
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++)
    h = max(h, P->h);
}

bool Solide::voisins_face(int num_face) {
  int part_1 = faces[num_face].voisins[0];
  int part_2 = faces[num_face].voisins[1];
  //cout << part_1 << " " << part_2 << endl;
  //cout << solide[part_1].impact_splitting << " " << solide[part_2].impact_splitting << endl;
  //cout << solide[part_1].from_splitting << " " << solide[part_2].from_splitting << endl;
  

  vector<int> tous_voisins; //Va Ãªtre rempli des voisins des 2 particules qui peuvent Ãªtre candidats pour former le tetra associÃ© Ã  la face !

  bool tetra_ok = false;
  if(not(faces[num_face].split)){
    
    // //Ancienne version : recherche sur les particules voisines des particules voisines de la face
    // for(std::vector<int>::iterator G=solide[part_1].faces.begin();G!=solide[part_1].faces.end();G++){
    //   //cout << "Num face teste : " << *G << endl;
    //   if(not(faces[*G] == faces[num_face]) && faces[*G].voisins[0] != part_1 && faces[*G].voisins[0] != -1)
    // 	tous_voisins.push_back(faces[*G].voisins[0]);
    //   else if(not(faces[*G] == faces[num_face]) && faces[*G].voisins[1] != part_1 && faces[*G].voisins[1] != -1)
    // 	tous_voisins.push_back(faces[*G].voisins[1]);
    // }
    // //Boucle sur deuxiÃ¨me groupe de particules
    // for(std::vector<int>::iterator G=solide[part_2].faces.begin();G!=solide[part_2].faces.end();G++){
    //   //cout << "Num face teste : " << *G << endl;
    //   if(not(faces[*G] == faces[num_face]) && faces[*G].voisins[0] != part_2 && faces[*G].voisins[0] != -1)
    // 	tous_voisins.push_back(faces[*G].voisins[0]);
    //   else if(not(faces[*G] == faces[num_face]) && faces[*G].voisins[1] != part_2 && faces[*G].voisins[1] != -1)
    // 	tous_voisins.push_back(faces[*G].voisins[1]);
    // }
    // //Tous_voisins rempli a priori
    // //cout << "Taille recherche : " << tous_voisins.size() << endl;

    //Nouvelle version : recherche sur les particules qui partagent un vertex avec la face
    for(std::vector<int>::iterator V=faces[num_face].vertex.begin();V!=faces[num_face].vertex.end();V++){
      for(std::vector<int>::iterator P=vertex[*V].particules.begin();P!=vertex[*V].particules.end();P++){
	if(*P!=part_1 && *P!=part_2){
	  tous_voisins.push_back(*P);
	  if(solide[*P].split){
	    cout << "particule splittee !!!" << endl;
	    getchar();
	  }
	}
      }
    }
    

    // //Test visualisation voisins
    // std::ofstream vtk;
    // vtk.open("voisins_face_test.vtk",ios::out);
    // //Initialisation du fichier vtk
    // vtk << "# vtk DataFile Version 3.0" << endl;
    // vtk << "#Test voisins_face" << endl;
    // vtk << "ASCII" << endl;
    // vtk<< "\n";
    // vtk << "DATASET UNSTRUCTURED_GRID" << endl;
    // vtk << "POINTS " << 4*(tous_voisins.size()+2) << " DOUBLE" << endl;
    // //Sortie des points
    // for(std::vector<int>::iterator p=tous_voisins.begin();p!=tous_voisins.end();p++){
    //   for(std::vector<int>::iterator V=(solide[*p].vertices).begin();V!=(solide[*p].vertices).end();V++){
    // 	vtk << solide[*p].mvt_t(vertex[*V].pos) << endl;
    //   }
    // }
    // for(std::vector<int>::iterator V=(solide[part_1].vertices).begin();V!=(solide[part_1].vertices).end();V++){
    // 	vtk << solide[part_1].mvt_t(vertex[*V].pos) << endl;
    // }
    // for(std::vector<int>::iterator V=(solide[part_2].vertices).begin();V!=(solide[part_2].vertices).end();V++){
    // 	vtk << solide[part_2].mvt_t(vertex[*V].pos) << endl;
    // }
    // //Sortie des cellules
    // vtk << "CELLS " << tous_voisins.size()+2 << " " << 5*(tous_voisins.size()+2) << endl;
    // int point_tmp = 0;
    // for(std::vector<int>::iterator p=tous_voisins.begin();p!=tous_voisins.end();p++){
    //   vtk << 4;
    //   for(std::vector<int>::iterator V=(solide[*p].vertices).begin();V!=(solide[*p].vertices).end();V++){
    // 	vtk << " " << point_tmp;
    // 	point_tmp++;
    //   }
    //   vtk << endl;
    // }
    // vtk << 4;
    // for(std::vector<int>::iterator V=(solide[part_1].vertices).begin();V!=(solide[part_1].vertices).end();V++){
    //   vtk << " " << point_tmp;
    //   point_tmp++;
    // }
    // vtk << endl;
    // vtk << 4;
    // for(std::vector<int>::iterator V=(solide[part_2].vertices).begin();V!=(solide[part_2].vertices).end();V++){
    //   vtk << " " << point_tmp;
    //   point_tmp++;
    // }
    // vtk << endl;
    // vtk << "CELL_TYPES " << tous_voisins.size()+2 << endl;
    // for(int p=0;p<(tous_voisins.size()+2);p++){
    //   vtk << 10 << endl;
    // }
    // vtk << "CELL_DATA " << tous_voisins.size()+2 << endl;
    // vtk << "SCALARS voisin double 1" << endl;
    // vtk << "LOOKUP_TABLE default" << endl;
    // for(int p=0;p<(tous_voisins.size());p++){
    //   vtk << 1 << endl;
    // }
    // vtk << 0 << endl;
    // vtk << 0 << endl;
    // //Fin visualisation test

    bool test_volume=false;
    double max_c_min=-10000.;
    double max_d_min=-10000.;
    for(std::vector<int>::iterator G=tous_voisins.begin();G!=tous_voisins.end() && not(tetra_ok);G++){
      for(std::vector<int>::iterator I=G + 1;I!=tous_voisins.end() && not(tetra_ok);I++){
	//for(std::vector<int>::iterator I=tous_voisins.begin();I!=tous_voisins.end() && not(tetra_ok);I++){
	int voisin1 = *G;
	int voisin2 = *I;
	double vol = cross_product(Vector_3(solide[part_1].x0,solide[part_2].x0),Vector_3(solide[part_1].x0,solide[voisin1].x0))*Vector_3(solide[part_1].x0,solide[voisin2].x0); //6 * volume du tetra associé à la face
	//cout << "Volume : " << vol << endl;
	//if(vol > pow(10., -20.)) { //Attention avec cette valeur !! Il faut savoir qu'elle est là ! //-10.
	if(std::abs(vol) > 0.){
	  test_volume = true;
	  double c_part_1 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_2].x0, solide[voisin2].x0), Vector_3(solide[part_2].x0, solide[voisin1].x0)) ) / vol;
	  double d_part_1 = c_part_1 * std::abs(vol) / sqrt(cross_product(Vector_3(solide[part_2].x0, solide[voisin2].x0), Vector_3(solide[part_2].x0, solide[voisin1].x0)).squared_length());
	  double c_part_2 = (Vector_3(solide[part_1].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0)) ) / vol;
	  double d_part_2 = c_part_2 * std::abs(vol) / sqrt(cross_product(Vector_3(solide[part_1].x0, solide[voisin1].x0), Vector_3(solide[part_1].x0, solide[voisin2].x0)).squared_length());
	  double c_voisin1 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_1].x0, solide[voisin2].x0), Vector_3(solide[part_1].x0, solide[part_2].x0)) ) / vol;
	  double d_voisin1 = c_voisin1 * std::abs(vol) / sqrt(cross_product(Vector_3(solide[part_1].x0, solide[voisin2].x0), Vector_3(solide[part_1].x0, solide[part_2].x0)).squared_length());
	  double c_voisin2 = (Vector_3(solide[part_2].x0, faces[num_face].centre) * cross_product(Vector_3(solide[part_1].x0, solide[part_2].x0), Vector_3(solide[part_1].x0, solide[voisin1].x0)) ) / vol;
	  double d_voisin2 = c_voisin2 * std::abs(vol) / sqrt(cross_product(Vector_3(solide[part_1].x0, solide[part_2].x0), Vector_3(solide[part_1].x0, solide[voisin1].x0)).squared_length());
	  
	  //cout << "Face " << num_face << " Coords bary : " << c_part_1 << " " << c_part_2 << " " << c_voisin1 << " " << c_voisin2 << endl;

	  // if(std::abs(c_part_1+c_part_2+c_voisin1+c_voisin2-1)>1e-10){
	  //   cout << "Face=" << num_face << " Volume : " << vol << endl;
	  //   cout << "Somme bary=" << c_part_1+c_part_2+c_voisin1+c_voisin2 << endl;
	  //   cout << "c_part_1=" << c_part_1 << " c_part_2=" << c_part_2 << " c_voisin1=" << c_voisin1 << " c_voisin2=" << c_voisin2 << endl;
	  //   //getchar();
	  // }
	  
	  double c_min=min(min(c_part_1,c_part_2),min(c_voisin1,c_voisin2));
	  double d_min=min(min(d_part_1,d_part_2),min(d_voisin1,d_voisin2));
	  
	  //cout << "c_min=" <<  c_min << endl;
	  //cout << "Face=" << num_face << " d_min=" << d_min << " c_min=" << c_min << endl;
	  //if(c_min>max_c_min){
	  if(d_min>max_d_min && (c_min>=0. || std::abs(vol)>1e-4*h)){
	    //cout << "chgt d_min=" << d_min << endl;
	    max_c_min = max(max_c_min,c_min);
	    max_d_min = max(max_d_min,d_min);
	    faces[num_face].c_reconstruction.clear();
	    faces[num_face].reconstruction.clear();
	    faces[num_face].c_reconstruction.push_back(c_part_1);
	    faces[num_face].c_reconstruction.push_back(c_part_2);
	    faces[num_face].c_reconstruction.push_back(c_voisin1);
	    faces[num_face].c_reconstruction.push_back(c_voisin2);
	    faces[num_face].reconstruction.push_back(part_1);
	    faces[num_face].reconstruction.push_back(part_2);
	    faces[num_face].reconstruction.push_back(voisin1);
	    faces[num_face].reconstruction.push_back(voisin2);
	    if(c_min>=0.){
	      tetra_ok = true;
	    }
	  }
	}
      }
    }

    // //Test visualisation voisins retenus dans la reconstruction
    // std::ofstream vtk_rec;
    // vtk_rec.open("reconstruction_test.vtk",ios::out);
    // //Initialisation du fichier vtk
    // vtk_rec << "# vtk DataFile Version 3.0" << endl;
    // vtk_rec << "#Test voisins_face" << endl;
    // vtk_rec << "ASCII" << endl;
    // vtk_rec<< "\n";
    // vtk_rec << "DATASET UNSTRUCTURED_GRID" << endl;
    // vtk_rec << "POINTS " << 16 << " DOUBLE" << endl;
    // //Sortie des points
    // for(std::vector<int>::iterator p=faces[num_face].reconstruction.begin();p!=faces[num_face].reconstruction.end();p++){
    //   for(std::vector<int>::iterator V=(solide[*p].vertices).begin();V!=(solide[*p].vertices).end();V++){
    // 	vtk_rec << solide[*p].mvt_t(vertex[*V].pos) << endl;
    //   }
    // }
    // //Sortie des cellules
    // vtk_rec << "CELLS " << faces[num_face].reconstruction.size() << " " << 5*(faces[num_face].reconstruction.size()) << endl;
    // point_tmp = 0;
    // for(std::vector<int>::iterator p=faces[num_face].reconstruction.begin();p!=faces[num_face].reconstruction.end();p++){
    //   vtk_rec << 4;
    //   for(std::vector<int>::iterator V=(solide[*p].vertices).begin();V!=(solide[*p].vertices).end();V++){
    // 	vtk_rec << " " << point_tmp;
    // 	point_tmp++;
    //   }
    //   vtk_rec << endl;
    // }
    // vtk_rec << "CELL_TYPES " << faces[num_face].reconstruction.size() << endl;
    // for(int p=0;p<faces[num_face].reconstruction.size();p++){
    //   vtk_rec << 10 << endl;
    // }
    // vtk_rec << "CELL_DATA " << faces[num_face].reconstruction.size() << endl;
    // vtk_rec << "SCALARS voisin double 1" << endl;
    // vtk_rec << "LOOKUP_TABLE default" << endl;
    // vtk_rec << 0 << endl << 0 << endl << 1 << endl << 1 << endl;
    // vtk_rec << "SCALARS coeff_reconstruction double 1" << endl;
    // vtk_rec << "LOOKUP_TABLE default" << endl;
    // for(std::vector<double>::iterator c=faces[num_face].c_reconstruction.begin();c!=faces[num_face].c_reconstruction.end();c++){
    //   vtk_rec << *c << endl;
    // }
    
    
    // //Fin visualisation test
    
    if(faces[num_face].c_reconstruction.size()!=4 || faces[num_face].reconstruction.size()!=4){
      cout << "Taille reconstruction=" << faces[num_face].c_reconstruction.size() << " " << faces[num_face].reconstruction.size() << endl;
      getchar();
    } else {
      double c_part_1 = faces[num_face].c_reconstruction[0];
      double c_part_2 = faces[num_face].c_reconstruction[1];
      double c_voisin1 = faces[num_face].c_reconstruction[2];
      double c_voisin2 = faces[num_face].c_reconstruction[3];
      int part_1 = faces[num_face].reconstruction[0];
      int part_2 = faces[num_face].reconstruction[1];
      int voisin1 = faces[num_face].reconstruction[2];
      int voisin2 = faces[num_face].reconstruction[3];
      if(std::abs(c_part_1+c_part_2+c_voisin1+c_voisin2-1)>1e-10){
	cout << "Face=" << num_face << " Somme bary=" << c_part_1+c_part_2+c_voisin1+c_voisin2 << endl;
	cout << "c_part_1=" << c_part_1 << " c_part_2=" << c_part_2 << " c_voisin1=" << c_voisin1 << " c_voisin2=" << c_voisin2 << endl;
	getchar();
      }
    }
    //cout << "max(c_min)=" << max_c_min << endl;
    if(max_c_min<0){
      cout << "max(c_min) negatif !" << endl;
      //getchar();
    }
    
    
    if(not(test_volume)){
      cout << "Tous les tetras sont aplatis !" << endl;
      getchar();
    }
    
      
  }
  
  return tetra_ok;
}

void Solide::Solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  /*for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++) {
    if(not(P->split))
      P->solve_position(dt, flag_2d, t, T);
  }*/
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++) {
    F->solve_position(dt, t, T);
  }
}

void Solide::Solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  //Mise a jour du coeff d'amortissement dans chaque particule
  //Force_damping(dt, Amort, t, T);
  //Predicteur
  /*for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(t < T/2.)
      P->solve_vitesse(dt, flag_2d, Amort, t , T);
    else
    P->solve_vitesse_predictor(dt, flag_2d, Amort, t , T);
    //P->solve_vitesse_MEMM(dt, flag_2d, Amort, t , T);
    P->solve_vitesse(dt, flag_2d, Amort, t , T);
  }*/
  //Calcul des forces d'amortissement
  //Force_damping(dt, Amort, t, T);
  //Correcteur
  /*for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    P->solve_vitesse_corrector(dt, flag_2d, Amort, t , T);
    }*/

  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++) {
    //F->solve_vitesse(dt, t, T);
    F->solve_vitesse_MEMM(dt, t , T);
  }
}

void Solide::test_fissuration() {
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++) {
    F->test_fissuration(Gc);
  }
}

void Solide::Forces(const int& N_dim, const double& dt, const double& t, const double& T){
  //Mise a zero de l'integrale des forces
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
    F->F = Vector_3(0.,0.,0.);
  }
  //Integration par points de Gauss
  //Point milieu
  double theta=0.5; //theta=1. pour intégration Verlet //Theta=0.5 pour MEMM
  double weight = 1.;
  Forces_internes_bis(dt, theta, weight, t, T);
}

void Solide::stresses_bis_MEMM(const double& theta, const double& t, const double& T){
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){ //Calcul des déformations et des contraintes
    P->discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    P->discrete_gradient.col2 = Vector_3(0., 0., 0.);
    P->discrete_gradient.col3 = Vector_3(0., 0., 0.);
    P->discrete_sym_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    P->discrete_sym_gradient.col2 = Vector_3(0., 0., 0.);
    P->discrete_sym_gradient.col3 = Vector_3(0., 0., 0.);
    for(int i=0 ; i < P->faces.size() ; i++){
      int f = P->faces[i];
      Vector_3 nIJ = faces[f].normale;
      if(faces[f].BC == 0 && nIJ * Vector_3(P->x0, faces[f].centre) < 0.)
	nIJ = -nIJ; //Normale pas dans le bon sens...
      Matrix Dij_n(tens(0.5 * (faces[f].I_Dx + faces[f].I_Dx_prev),  nIJ));
      Matrix Dij_n_sym(tens_sym(0.5 * (faces[f].I_Dx + faces[f].I_Dx_prev),  nIJ));
      P->discrete_gradient += faces[f].S /  P->V * Dij_n;
      P->discrete_sym_gradient += faces[f].S /  P->V * Dij_n_sym;
    }

    //calcul des contraintes
    P->contrainte = lambda * (P->discrete_sym_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_sym_gradient - P->epsilon_p);

    //Plastification si on dépasse le critère  
    P->seuil_elas = A; // + B * pow(P->def_plas_cumulee, n);
    if((P->contrainte - H * P->epsilon_p).VM() > P->seuil_elas) { //On sort du domaine élastique.
      Matrix n_elas( 1. / ((P->contrainte).dev()).norme() * (P->contrainte).dev() ); //Normale au domaine élastique de Von Mises
      double delta_p = ((P->contrainte - H * P->epsilon_p).VM() - A) / (2*mu + H);
      //P->def_plas_cumulee += delta_p;
      //cout << "delta_p : " << delta_p << endl;
      //P->epsilon_p += delta_p * n_elas;
      //cout << "Def plas : " << P->epsilon_p.col1 << endl;
      //P->contrainte = lambda * (P->discrete_sym_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_sym_gradient - P->epsilon_p); //Recalcul des contraintes après plastification
    }
  }
}

void Solide::stresses_bis(const double& theta, const double& t, const double& T){
  //for(std::vector<Face>::iterator F=faces.begin(); F!=faces.end(); F++) { //On impose la velur des DDL sur des faces Neumann
    /*if(F->BC == 1) { //Dirichlet
      F->I_Dx = displacement_BC_bis(F->centre, Vector_3(0.,0.,0.), t, T) * F->normale;
    }
    else if(F->BC == -1) //On impose la contrainte voulue sur les faces Neumann
      F->F = F->F + //Contrainte Neumann sur la face !!! */

    /*double def_ref = 0.001; // * t / T; //Solution statique
    //F->I_Dx = def_ref * F->centre.z() * F->normale;
    F->I_Dx.vec[2] = F->centre.z() * def_ref;
    F->I_Dx.vec[0] = -0.3 * F->centre.x() * def_ref;
    F->I_Dx.vec[1] = -0.3 * F->centre.y() * def_ref;*/
  //}

  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){ //Calcul des déformations et des contraintes
    P->discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    P->discrete_gradient.col2 = Vector_3(0., 0., 0.);
    P->discrete_gradient.col3 = Vector_3(0., 0., 0.);
    P->discrete_sym_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    P->discrete_sym_gradient.col2 = Vector_3(0., 0., 0.);
    P->discrete_sym_gradient.col3 = Vector_3(0., 0., 0.);
    for(int i=0 ; i < P->faces.size() ; i++){
      int f = P->faces[i];
      Vector_3 nIJ = faces[f].normale;
      if(faces[f].BC == 0 && nIJ * Vector_3(P->x0, faces[f].centre) < 0.)
	nIJ = -nIJ; //Normale pas dans le bon sens...
      if(not(faces[f].fissure)) {
	Matrix Dij_n(tens(faces[f].I_Dx,  nIJ));
	Matrix Dij_n_sym(tens_sym(faces[f].I_Dx,  nIJ));
	P->discrete_gradient += faces[f].S /  P->V * Dij_n;
	P->discrete_sym_gradient += faces[f].S /  P->V * Dij_n_sym;
      }
      else if(faces[f].fissure) {
	int voi;
	if(P->id == faces[f].voisins[0])
	  voi = 0;
	else if(P->id == faces[f].voisins[1])
	  voi = 1;
	Matrix Dij_n(tens(faces[f].Dx[voi],  nIJ));
	Matrix Dij_n_sym(tens_sym(faces[f].Dx[voi],  nIJ));
	P->discrete_gradient += faces[f].S /  P->V * Dij_n;
	P->discrete_sym_gradient += faces[f].S /  P->V * Dij_n_sym;
      }
    }

    //calcul des contraintes
    P->contrainte = lambda * (P->discrete_sym_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_sym_gradient - P->epsilon_p);

    //Plastification si on dépasse le critère  
    P->seuil_elas = A; // + B * pow(P->def_plas_cumulee, n);
    if((P->contrainte - H * P->epsilon_p).VM() > P->seuil_elas) { //On sort du domaine élastique.
      Matrix n_elas( 1. / ((P->contrainte).dev()).norme() * (P->contrainte).dev() ); //Normale au domaine élastique de Von Mises
      double delta_p = ((P->contrainte - H * P->epsilon_p).VM() - A) / (2*mu + H);
      //P->def_plas_cumulee += delta_p;
      //cout << "delta_p : " << delta_p << endl;
      //P->epsilon_p += delta_p * n_elas;
      //cout << "Def plas : " << P->epsilon_p.col1 << endl;
      //P->contrainte = lambda * (P->discrete_sym_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_sym_gradient - P->epsilon_p); //Recalcul des contraintes après plastification
    }
  }
}

void Solide::stresses(const double& theta, const double& t, const double& T){ //Calcul de la contrainte dans toutes les particules ; theta indique qu'on calcule la force a partir de la position au temps t+theta*dt
  for(int i=0; i<faces.size(); i++){ //Calcul de la reconstruction sur chaque face
    if(faces[i].BC == 0) //cad face dans bulk et donc I_Dx reconstruit
      faces[i].I_Dx = Vector_3(0., 0., 0.); //Remise Ã  zÃ©ro. Si particule sur le bord, on a bien I_Dx = (0., 0., 0.)
    //cout << "BC : " << faces[i].BC << endl;
    //Vector_3 test_pos(0., 0., 0.);
    if(faces[i].BC == 1) { //Dirichlet
      faces[i].I_Dx = Vector_3(0., 0., 0.);
      double def_ref = 0.001; // * t / T; //Solution statique
      faces[i].I_Dx = def_ref * faces[i].centre.z() * faces[i].normale;
      //faces[i].I_Dx = faces[i].I_Dx + displacement_BC_bis(faces[i].centre, solide[faces[i].voisins[0]].Dx, t, T) * faces[i].normale;
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
    else if(faces[i].BC == 0 && not(faces[i].split)) { //Cad particule dans le bulk. Donc reconstruction !
      for(int j=0; j<faces[i].reconstruction.size() ; j++) {
	Vector_3 Dx = theta*solide[faces[i].reconstruction[j]].Dx + (1.-theta)*solide[faces[i].reconstruction[j]].Dxprev;
	faces[i].I_Dx = faces[i].I_Dx + faces[i].c_reconstruction[j] * Dx;
	//cout << "Interpolation : " << faces[i].I_Dx << endl;
      }
      //faces[i].I_Dx = (solide[faces[i].voisins[0]].Dx + solide[faces[i].voisins[1]].Dx) / 2.;

      /*if(abs(faces[i].I_Dx[2] - 4. / 3. * faces[i].centre.z()) > 0.00001)
	cout << "ProblÃ¨me reconstruction sur face : " << i << endl;*/
      /*if(sqrt((test_pos - Vector_3(Point_3(0.,0.,0.),faces[i].centre)).squared_length()) > pow(10., -11.))
	cout << "ProblÃ¨me reconstruction barycentre face : " << i << endl;*/
    }
  }
  
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      bool test_continuer;
      int nb_iterations = 0;
      do {
	test_continuer = false;
	nb_iterations++;
	P->discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
	P->discrete_gradient.col2 = Vector_3(0., 0., 0.);
	P->discrete_gradient.col3 = Vector_3(0., 0., 0.);
	P->discrete_sym_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
	P->discrete_sym_gradient.col2 = Vector_3(0., 0., 0.);
	P->discrete_sym_gradient.col3 = Vector_3(0., 0., 0.);
	std::vector<int> num_faces; //Sert à récupérer le numéro des faces avec BC de Neumann si nécessaire
	for(int i=0 ; i < P->faces.size() ; i++){
	  int f = P->faces[i];
	  if(faces[f].BC != 0) { //car on ne fait ces calculs sur toutes les faces au bord
	    num_faces.push_back(f);
	  }
	  Vector_3 nIJ = faces[f].normale;
	  if(faces[f].BC == 0 && nIJ * Vector_3(P->x0, faces[f].centre) < 0.)
	    nIJ = -nIJ; //Normale pas dans le bon sens...
	  Matrix Dij_n(tens(faces[f].I_Dx,  nIJ));
	  Matrix Dij_n_sym(tens_sym(faces[f].I_Dx,  nIJ));
	  if(faces[f].BC >= 0){ //On ajoute pas les faces de Neumann ni de Dirichlet car on doit recalculer la valeur sur la face// >= 0 avant test
	    P->discrete_gradient += faces[f].S /  P->V * Dij_n;
	    P->discrete_sym_gradient += faces[f].S /  P->V * Dij_n_sym;
	  }
	}

	P->contrainte = lambda * (P->discrete_sym_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_sym_gradient - P->epsilon_p); //Premier calcul pour calcul des déplacements sur bords de Neumann
      
	//On reconstruit la valeur du déplacement sur les faces de Neumann
	if(num_faces.size() > 0) {
	  //cout << "Num particule : " << P->id << endl;
	  /*cout << "Apres calcul : " << endl;
	    cout << P->contrainte.col1 << endl;
	  cout << P->contrainte.col2 << endl;
	  cout << P->contrainte.col3 << endl;*/
	  //reconstruction_faces_neumann(num_faces, P->contrainte, t, P->V, T); //Calcul la valeur des déplacements sur faces de Neumann
	  
	  //On recalcul le gradient discret
	  P->discrete_sym_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
	  P->discrete_sym_gradient.col2 = Vector_3(0., 0., 0.);
	  P->discrete_sym_gradient.col3 = Vector_3(0., 0., 0.);
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
	    P->discrete_sym_gradient += faces[f].S /  P->V * Dij_n;
	    P->discrete_gradient += faces[f].S / P->V * tens(faces[f].I_Dx, nIJ);
	  }
	}
	
	P->contrainte = lambda * (P->discrete_sym_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_sym_gradient - P->epsilon_p); //Contrainte complète
	
	/*for(int i=0 ; i < P->faces.size() ; i++){ //On recalcule les contraintes avec toutes les contributions
	  int f = P->faces[i];
	  if(faces[f].BC == -1 && sqrt((P->contrainte * faces[f].normale).squared_length()) > 1.)
	  cout << "ProblÃ¨me Neumann homogÃ¨ne : " << sqrt((P->contrainte * faces[f].normale).squared_length()) << ". Num face : " << faces[f].id << endl;
	  }*/
	
	//Plastification si on dépasse le critère  
	P->seuil_elas = A; // + B * pow(P->def_plas_cumulee, n);
	if((P->contrainte - H * P->epsilon_p).VM() > P->seuil_elas) { //On sort du domaine élastique.
	  //while( (P->contrainte).VM() > A) { //On dépasse le critère plastique on refait un return mapping pour essayer de converger
	  //cout << "Plastification !!!!" << endl;
	  //Plastification
	  Matrix n_elas( 1. / ((P->contrainte).dev()).norme() * (P->contrainte).dev() ); //Normale au domaine élastique de Von Mises
	  double delta_p = ((P->contrainte - H * P->epsilon_p).VM() - A) / (2*mu + H);
	  P->def_plas_cumulee += delta_p;
	  //cout << "delta_p : " << delta_p << endl;
	  P->epsilon_p += delta_p * n_elas;
	  //cout << "Def plas : " << P->epsilon_p.col1 << endl;
	  P->contrainte = lambda * (P->discrete_sym_gradient - P->epsilon_p).tr() * unit() + 2*mu * (P->discrete_sym_gradient - P->epsilon_p); //Recalcul des contraintes après plastification
	}
	if(num_faces.size() > 0) { //On vérifie qu'on a toujours les bonnes BC de Neumann
	  for(int i=0 ; i < num_faces.size() ; i++) {
	    if( faces[i].BC == -1 && sqrt((P->contrainte * faces[i].normale).squared_length()) > 1.) {
	      test_continuer = true;
	      //cout << "On continue la boucle 1 !" << endl;
	      break;
	    }
	    //else if( sqrt( (P->contrainte * faces[i].normale - ((P->contrainte * faces[i].normale) * faces[i].normale) * faces[i].normale).squared_length()) > 1. && faces[i].BC == 1) {
	    else if( faces[i].BC == 1 && sqrt(pow((P->contrainte * faces[i].normale) * faces[i].vec_tangent_1, 2.)) > 1. && sqrt(pow((P->contrainte * faces[i].normale) * faces[i].vec_tangent_2, 2.)) > 1. ) {
	      test_continuer = true;
	      /*cout << "On continue la boucle 2 !" << endl;
	      cout << sqrt(pow((P->contrainte * faces[i].normale) * faces[i].vec_tangent_1, 2.)) << " " << sqrt(pow((P->contrainte * faces[i].normale) * faces[i].vec_tangent_2, 2.)) << endl;
	      cout << P->epsilon_p.col1 << endl;
	      cout << P->epsilon_p.col2 << endl;
	      cout << P->epsilon_p.col3 << endl;
	      cout << P->def_plas_cumulee << endl;*/
	      break;
	    }
	  }
	}
      } while((P->contrainte - H * P->epsilon_p).VM() > P->seuil_elas && test_continuer && nb_iterations < 10); //nb_iterations < 10
      //if(nb_iterations == 10)
      //cout << "Particule : " << P->id << " pas de convergence entre plasticite et BC !" << endl;
      /*if((P->contrainte).VM() > P->seuil_elas)
	cout << "Von Mises : " << (P->contrainte - H * P->epsilon_p).VM() << endl;*/
    }
  }
}

void Solide::Forces_internes_bis(const double& dt, const double& theta, const double& weight, const double& t, const double& T){ //Calcul des forces pour chaque particule ; theta indique qu'on calcule la force a partir de la position au temps t+theta*dt
  if(theta > 0.9) //Cad integration Verlet
    stresses_bis(theta, t, T);
  else if(theta > 0.4 && theta < 0.6) //Integration MEMM
    stresses_bis_MEMM(theta, t, T);

  for(std::vector<Face>::iterator F=faces.begin(); F!=faces.end(); F++) { //Vraies forces
    int voisin1 = F->voisins[0];
    int voisin2 = F->voisins[1];
    F->F = -F->S * solide[voisin1].contrainte * F->normale; //-
    F->Forces[0] = -F->S * solide[voisin1].contrainte * F->normale;
    if(voisin2 >= 0) { //cad face pas sur un bord
      F->Forces[1] = F->S * solide[voisin2].contrainte * F->normale;
      F->F = F->F + F->S * solide[voisin2].contrainte * F->normale; //+
    }
  }
}


void Solide::Forces_internes(const double& dt, const double& theta, const double& weight, const double& t, const double& T){ //Calcul des forces pour chaque particule ; theta indique qu'on calcule la force a partir de la position au temps t+theta*dt
  stresses(theta, t, T);
  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++) //Remet Ã  zÃ©ro toutes les forces
    P->Fi = Vector_3(0.,0.,0.);
  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++){
    if(not(P->split)) {
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
	  double signe = 0.;
	  if(nIJ * Vector_3(P->x0, faces[num_face].centre) < 0.) {
	    nIJ = -nIJ; //Normale pas dans le bon sens...
	    signe = -1.;
	  }
	  else
	    signe = 1.;

	  //P->Fi = P->Fi +  2 * mu * (solide[voisin].Dx - P->Dx); //OK flux à 2 points
	  //P->Fi = P->Fi + faces[num_face].S * (solide[part_1].contrainte + solide[part_2].contrainte) / 2. * nIJ; //OK Voronoi
	  P->Fi = P->Fi + faces[num_face].S * P->contrainte * nIJ; //Force sur la particule sommée sur toutes les faces. Pour Tetra !
	  solide[aux_1].Fi = solide[aux_1].Fi - faces[num_face].S * c_aux_1 * P->contrainte * nIJ; //Force sur autres particules en interaction avec P
	  solide[aux_2].Fi = solide[aux_2].Fi - faces[num_face].S * c_aux_2 * P->contrainte * nIJ;
	  solide[aux_3].Fi = solide[aux_3].Fi - faces[num_face].S * c_aux_3 * P->contrainte * nIJ;
	  solide[aux_4].Fi = solide[aux_4].Fi - faces[num_face].S * c_aux_4 * P->contrainte * nIJ;

	  //Ajout de la force issue de la pénalisation
	  P->Fi = P->Fi + faces[num_face].S * eta / faces[num_face].h * (solide[faces[num_face].voisins[0]].u - solide[faces[num_face].voisins[1]].u) * signe; //signe pour changer le sens du saut selon la particule choisie
	}
	else if(faces[num_face].BC == -1) { //Calcul forces sur DDL bords Dirichlet. Vaut 0 exactement en Neumann Homogène. A vérifier... // == 1
	  //cout << "Face au bord" << endl;
	  int part = faces[num_face].voisins[0];
	  Vector_3 nIJ = faces[num_face].normale;
	  /*if((faces[num_face].S * solide[part].contrainte * nIJ).squared_length() > 1.)
	    cout << "Pb avec force sur bord de Neumann : " << (faces[num_face].S * solide[part].contrainte * nIJ).squared_length() << endl;*/
	  P->Fi = P->Fi + faces[num_face].S * solide[part].contrainte * nIJ; //pow(10., 7.) * nIJ;
	  //cout << "Force volume : " << P->Fi * nIJ << endl;
	}
	else if(faces[num_face].BC == 1) { //Calcul forces sur DDL bords Dirichlet. Vaut 0 exactement en Neumann Homogène. A vérifier... // == 1
	  //cout << "Face au bord" << endl;
	  int part = faces[num_face].voisins[0];
	  Vector_3 nIJ = faces[num_face].normale;
	  P->Fi = P->Fi + faces[num_face].S * ((solide[part].contrainte * nIJ) * nIJ) * nIJ; //On ne met que le composante imposée
	}
      }
    }
  }

  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++) //Pour debug de la solution statique
    cout << "Particule :" << P->id << " Force : " << P->Fi << endl;
  
  //Mise a jour de l'integrale en temps des forces
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    P->Fi_int = P->Fi_int + weight*P->Fi;
  }
}

void Solide::Force_damping(const double& dt, const double& Amort, const double& t, const double& T)
{//Calcul des forces d'amortissement pour chaque particule
  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++){ //Remet à zéro des forces d'amortissement et de coefficients d'amortissement par lien
    P->F_damp = Vector_3(0.,0.,0.);
    P->damping = 0.;
  }
  
  for(std::vector<Particule>::iterator P=solide.begin(); P!=solide.end(); P++){
    for(int i=0 ; i<P->faces.size() ; i++){
      int num_face = P->faces[i]; //NumÃ©ro de la face dans l'ensemble des faces contenu dans le solide
      /*int part_1 = faces[num_face].voisins[0];
	int part_2 = faces[num_face].voisins[1];*/
      if(faces[num_face].BC == 0) { //Forces sur faces internes
	int part = faces[num_face].voisins[0];
	int part2 = faces[num_face].voisins[1];
	if(part==P->id){
	  //cout << "part=" << part << " P.id=" << P->id << endl;
	  part = faces[num_face].voisins[1];
	  part2 = faces[num_face].voisins[0];
	}
	const Particule& P1 = solide[part];
	if(part2!=P->id || part!=P1.id){
	  cout << "part2=" << part2 << " P.id=" << P->id << " part=" << part << " P1.id=" << P1.id << endl;
	  getchar();
	}
	double damping_link = Amort*sqrt(faces[num_face].S*(lambda+2*mu)*(P->m/P->V+P1.m/P1.V));
	P->F_damp = P->F_damp + damping_link * P1.u;
	P->damping += damping_link;
      }
    }
  }
}


const double Solide::Energie(){
  return Energie_cinetique()+Energie_potentielle();
}

const double Solide::Energie_MEMM(const double &t, const double &T){
  return Energie_cinetique_MEMM()+Energie_potentielle_MEMM(t, T);
}

const double Solide::Energie_cinetique(){ //Energie pas adaptée à MEMM non ?
  double E = 0.;
  /*for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split))
      E += 1./2. * P->m * (P->u + P->u_prev) * (P->u + P->u_prev) / 4.;
    //E += 1./2. * P->m * P->u * P->u;
    if(1./2. * P->m * (P->u + P->u_prev) * (P->u + P->u_prev) / 4. < 0.) {
      //cout << "Particule : " << P->id << " energie cinétique negative" << endl;
      //cout << "Volume : " << P->V << endl;
    }
    }*/

  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
    E += 0.5 * F->m * (F->u + F->u_prev) * (F->u + F->u_prev) / 4.;
  }
  
  return E;
}

const double Solide::Energie_cinetique_MEMM(){ //Energie cinétique adaptée à MEMM
  double E = 0.;
  for(std::vector<Face>::iterator F=faces.begin();F!=faces.end();F++){
    E += 0.5 * F->m * (F->u * F->u_prev);
  }
  
  return E;
}

const double Solide::Energie_potentielle(){
  double Ep = 0.;

  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    Ep += 0.5 * contraction_double(P->contrainte, P->discrete_sym_gradient - P->epsilon_p) * P->V /*+ B * pow(((P->second)).def_plas_cumulee, 1. + n) / (n + 1.)*/ + A * P->def_plas_cumulee * P->V + 0.5 * contraction_double(H * P->epsilon_p, P->epsilon_p) * P->V;
  }

  return Ep;
}

const double Solide::Energie_potentielle_MEMM(const double &t, const double &T){ //On recalcule les contraintes et def au pas de temps entier pour avoir la conservation de l'énergie

  stresses_bis(1., t, T); //Recalcul des contraintes au pas de temps entier

  double Ep = 0.;

  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    Ep += 0.5 * contraction_double(P->contrainte, P->discrete_sym_gradient - P->epsilon_p) * P->V /*+ B * pow(((P->second)).def_plas_cumulee, 1. + n) / (n + 1.)*/ + A * P->def_plas_cumulee * P->V + 0.5 * contraction_double(H * P->epsilon_p, P->epsilon_p) * P->V;
    //}
  }

  return Ep;
}

double Solide::pas_temps(const double& t, const double& T, const double& cfls, const double& E, const double& nu, const double& rhos){
  double eps = 1e-14;//std::numeric_limits<double>::epsilon();
  double dt = std::numeric_limits<double>::infinity();
  //Restriction CFL liee aux forces internes
  double cs = sqrt(E*(1.-nu)/rhos/(1.+nu)/(1.-2.*nu));
  //Calcul du rayon de la sphÃ¨re inscrite
  double rmin = std::numeric_limits<double>::infinity();
  for(int i=0;i<size();i++){
    Particule& P = solide[i];
    if(not(P.split)){
      P.calcul_sphere_inscrite(this);
      rmin = min(rmin,P.r);
    }
  }
  
  // double sigma = 100000.;
  // for(int i=0;i<size();i++){
  //   for(int j=0;j<solide[i].faces.size();j++){
  //     sigma = min(sigma,faces[solide[i].faces[j]].D0);
  //   }
  // }
  // cout << "sigma=" << sigma << endl;
  // for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
  //   for(int j=0;j<P->faces.size();j++){
  //     if(faces[P->faces[j]].voisins[0] >=0 && faces[P->faces[j]].voisins[1] >= 0){
  // 	dt = min(dt,cfls*faces[P->faces[j]].D0/cs);
  //     }
  //   }
  // }
  //cout << "cs=" << cs << endl;
  //getchar();
  dt = cfls*rmin/cs;
  //cout << "dt prescrit=" << dt << endl;
  //getchar();
  
  dt = min(dt,T-t);
  return dt;
}

void Solide::Impression(const int &n){ //Sortie au format vtk
  int nb_points = 0; //solide.size();
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){ //On ne compte pas les particules splitées
    if(not(P->split)){
      //nb_points += P->faces.size() * 3;
      //Si tetra
      if(P->vertices.size()==4){
	nb_points += P->vertices.size();
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  nb_points += 3;
	}
      }
    }
  }
  int nb_faces = 0;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){ //On ne compte pas les particules splitées
    if(not(P->split)){
      nb_faces += P->faces.size();
    }
  }
  int nb_part = 0;
  int size = 0;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){ //On ne compte pas les particules splitées
    if(not(P->split)){
      //Si tetra
      if(P->vertices.size()==4){
	nb_part += 1;
	size += 5;
      }  else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  nb_part += 1;
	  size += 4;
	}
      }     
    }
  }
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
  //int nb_points = 3 * 4 * nb_part;
  //int nb_points = 4*nb_part;
  //int nb_faces = 4 * nb_part;
  //int size = 3 * nb_faces;
  /*for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      size += F->nb_vertex; //Egale Ã  nb_points au final ?
  }*/
  //size += nb_faces; //Pk ? MystÃ¨re...
    
  //Initialisation du fichier vtk
  vtk << "# vtk DataFile Version 3.0" << endl;
  vtk << "#Simulation Euler" << endl;
  vtk << "ASCII" << endl;
  vtk<< "\n";
  vtk << "DATASET UNSTRUCTURED_GRID" << endl;
  vtk << "POINTS " << nb_points << " DOUBLE" << endl;

  //Sortie des points
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      //for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
      //for(std::vector<int>::iterator V=faces[*F].vertex.begin();V!=faces[*F].vertex.end();V++)
      //Si tetra
      if(P->vertices.size()==4){
	for(std::vector<int>::iterator V=(P->vertices).begin();V!=(P->vertices).end();V++){
	  vtk << P->mvt_t(vertex[*V].pos) << endl;
	}
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  for(std::vector<int>::iterator V=faces[*F].vertex.begin();V!=faces[*F].vertex.end();V++){
	    vtk << P->mvt_t(vertex[*V].pos) << endl;
	  }
	}
      }     
    }
  }
  vtk << endl;
    
  // //Sortie des faces
  // int point_tmp=0;
  // vtk << "CELLS " << nb_faces << " " << size << endl;
  // int compteur_vertex = 0;
  // for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
  //   for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
  //     vtk << faces[*F].vertex.size();
  //     for(int k=0 ; k<faces[*F].vertex.size() ; k++)
  // 	vtk << " " << compteur_vertex + k; //(faces[*F].vertex)[k];
  //     vtk << endl;
  //     compteur_vertex += faces[*F].vertex.size();
  //   }
  //   //vtk << endl;
  // }
  // vtk << endl;
  //Sortie des particules
  vtk << "CELLS " << nb_part << " " << size << endl;
  int point_tmp=0;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      //Si tetra
      if(P->vertices.size()==4){
	vtk << 4;
	for(std::vector<int>::iterator V=P->vertices.begin();V!=P->vertices.end();V++){
	  vtk << " " << point_tmp;
	  point_tmp++;
	}
	vtk << endl;
      }
      else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << 3;
	  for(std::vector<int>::iterator V=faces[*F].vertex.begin();V!=faces[*F].vertex.end();V++){
	    vtk << " " << point_tmp;
	    point_tmp++;
	  }
	  vtk << endl;
	}
      }  
    }
  }
  // vtk << "CELL_TYPES " << nb_faces << endl;
  // for(int i=0;i<nb_faces;i++){
  //   //vtk << 9 << endl; //Pour Quad
  //   vtk << 5 << endl; //Pour triangle
  // }
  vtk << "CELL_TYPES " << nb_part << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      //Si tetra
      if(P->vertices.size()==4){
	vtk << 10 << endl; //Pour le tÃ©tra
      }
      else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << 5 << endl;
	}
      }
    }
  }
  vtk << "\n";
  vtk << "POINT_DATA " << nb_points << endl;
  //Deplacement reconstruit
  vtk << "VECTORS deplacement_rec double" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      //Si tetra
      if(P->vertices.size()==4){
	for(std::vector<int>::iterator V=P->vertices.begin();V!=P->vertices.end();V++){
	  vtk << faces[P->faces[0]].I_Dx + P->discrete_gradient*(vertex[*V].pos - faces[P->faces[0]].centre) << endl;
	}
      }
      else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  for(std::vector<int>::iterator V=faces[*F].vertex.begin();V!=faces[*F].vertex.end();V++){
	    vtk << P->Dx + P->discrete_gradient*(vertex[*V].pos-P->x0) << endl;
	  }
	}
      }   
    }
  }
  vtk << "\n";
  //vtk << "CELL_DATA " << nb_faces << endl;
  vtk << "CELL_DATA " << nb_part << endl;
  //Deplacement
  vtk << "VECTORS deplacement double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    //for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
    if(not(P->split)) {
      //Si tetra
      if(P->vertices.size()==4){
	//vtk << P->Dx << endl;
	Vector_3 dep = Vector_3(0.,0.,0.);
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  dep = dep + faces[*F].I_Dx / 4.;
	}
	vtk << dep << endl;
      }
      else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << P->Dx << endl;
	}
      } 
    }
  }
  vtk << "\n";
  //Vitesse
  vtk << "VECTORS vitesse double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      //for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      //Si tetra
      if(P->vertices.size()==4){
	Vector_3 vit = Vector_3(0.,0.,0.);
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vit = vit + faces[*F].u / 4.;
	}
	vtk << vit << endl;
	//vtk << P->u << endl;
      }
      else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << P->u << endl;
	}
      }     
    }
  }
  vtk << "\n";
  /*Normale
  // vtk << "VECTORS normale double" << endl;
  // //vtk << "LOOKUP_TABLE default" << endl;
  // for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
  //   for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
  //   vtk << faces[*F].normale << endl;
  // }
  // vtk << "\n";
  //Tangent 1
  vtk << "VECTORS tan_1 double" << endl;
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
    if(not(P->split)) {
//for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
      //Si tetra
      if(P->vertices.size()==4){
	vtk << P->contrainte.col1[0] << " " << P->contrainte.col1[1] << " " << P->contrainte.col1[2] << endl;
	vtk << P->contrainte.col2[0] << " " << P->contrainte.col2[1] << " " << P->contrainte.col2[2] << endl;
	vtk << P->contrainte.col3[0] << " " << P->contrainte.col3[1] << " " << P->contrainte.col3[2] << endl;
      }
      else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << P->contrainte.col1[0] << " " << P->contrainte.col1[1] << " " << P->contrainte.col1[2] << endl;
	vtk << P->contrainte.col2[0] << " " << P->contrainte.col2[1] << " " << P->contrainte.col2[2] << endl;
	vtk << P->contrainte.col3[0] << " " << P->contrainte.col3[1] << " " << P->contrainte.col3[2] << endl;
	}
      }
    }
  }
  vtk << "\n";
  //DÃ©formations
  vtk << "TENSORS deformations double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      //Si tetra
      if(P->vertices.size()==4){
	//for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
	vtk << P->discrete_sym_gradient.col1[0] << " " << P->discrete_sym_gradient.col1[1] << " " << P->discrete_sym_gradient.col1[2] << endl;
	vtk << P->discrete_sym_gradient.col2[0] << " " << P->discrete_sym_gradient.col2[1] << " " << P->discrete_sym_gradient.col2[2] << endl;
	vtk << P->discrete_sym_gradient.col3[0] << " " << P->discrete_sym_gradient.col3[1] << " " << P->discrete_sym_gradient.col3[2] << endl;
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << P->discrete_sym_gradient.col1[0] << " " << P->discrete_sym_gradient.col1[1] << " " << P->discrete_sym_gradient.col1[2] << endl;
	vtk << P->discrete_sym_gradient.col2[0] << " " << P->discrete_sym_gradient.col2[1] << " " << P->discrete_sym_gradient.col2[2] << endl;
	vtk << P->discrete_sym_gradient.col3[0] << " " << P->discrete_sym_gradient.col3[1] << " " << P->discrete_sym_gradient.col3[2] << endl;
	}
      }
    }
  }
  vtk << "\n";
  //Epsilon_p
  vtk << "TENSORS epsilon_p double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      //for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
      //Si tetra
      if(P->vertices.size()==4){
	vtk << P->epsilon_p.col1[0] << " " << P->epsilon_p.col1[1] << " " << P->epsilon_p.col1[2] << endl;
	vtk << P->epsilon_p.col2[0] << " " << P->epsilon_p.col2[1] << " " << P->epsilon_p.col2[2] << endl;
	vtk << P->epsilon_p.col3[0] << " " << P->epsilon_p.col3[1] << " " << P->epsilon_p.col3[2] << endl;
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << P->epsilon_p.col1[0] << " " << P->epsilon_p.col1[1] << " " << P->epsilon_p.col1[2] << endl;
	vtk << P->epsilon_p.col2[0] << " " << P->epsilon_p.col2[1] << " " << P->epsilon_p.col2[2] << endl;
	vtk << P->epsilon_p.col3[0] << " " << P->epsilon_p.col3[1] << " " << P->epsilon_p.col3[2] << endl;
	}
      }
    }
  }
  vtk << "\n";
  //Deformation plastique cumulÃ©e
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
    if(not(P->split)) {
      //Si tetra
      if(P->vertices.size()==4){
	//for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
	vtk << (P->contrainte).VM() << endl;
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << (P->contrainte).VM() << endl;
	}
      }
    }
  }
  vtk << "\n";
  //NumÃ©ro des particules
  /*vtk << "SCALARS id double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++)
      vtk << P->id << endl;
  }
  vtk << endl;*/
  //Nombre de faces avec des BCs
  vtk << "SCALARS nb_faces_BC double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      std::vector<int> num_faces;
      for(int i=0 ; i < P->faces.size() ; i++){
	int f = P->faces[i];
	if(faces[f].BC != 0) { //car on ne fait ces calculs sur toutes les faces au bord
	  num_faces.push_back(f);
	}
      }
      //Si tetra
      if(P->vertices.size()==4){
	vtk << num_faces.size() << endl;
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << num_faces.size() << endl;
	}
      }
    }
  }
  //Particules issues du splitting
  vtk << "SCALARS from_splitting double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      //Si tetra
      if(P->vertices.size()==4){
	vtk << P->from_splitting << endl;
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << P->from_splitting << endl;
	}
      }
    }
  }
  //Particules impactees par le splitting
  vtk << "SCALARS impact_splitting double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      //Si tetra
      if(P->vertices.size()==4){
	vtk << P->impact_splitting << endl;
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << P->impact_splitting << endl;
	}
      }
    }
  }
  vtk << "SCALARS rayon_sphere double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)) {
      //Si tetra
      if(P->vertices.size()==4){
	vtk << P->r << endl;
      } else {
	//Particule polyedrale : sortie des faces triangulaires
	for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	  vtk << P->r << endl;
	}
      }
    }
  }
  vtk.close();
}

void Solide::Impression_faces(const int &n){ //Sortie au format vtk des interpolations aux faces
  int nb_part = 0;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){ //On ne compte pas les particules splitées
    if(not(P->split)){
      nb_part += 1;
    }
  }
  std::ostringstream oss;
  oss << "faces" << n << ".vtk";
  string s = oss.str();
  const char* const solidevtk = s.c_str();
    
  //Ouverture des flux en donne en ecriture
  std::ofstream vtk;
  vtk.open(solidevtk,ios::out);
  if(not(vtk.is_open()))
    cout <<"ouverture de faces" << n << ".vtk rate" << endl;
  //vtk << setprecision(15);
  
  //Pour tetras !
  int nb_faces = 0;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){ //On ne compte pas les particules splitées
    if(not(P->split)){
      nb_faces += P->faces.size();
    }
  }
  int nb_points = 3 * nb_faces;
  int size = 4 * nb_faces;
    
  //Initialisation du fichier vtk
  vtk << "# vtk DataFile Version 3.0" << endl;
  vtk << "#Simulation Euler" << endl;
  vtk << "ASCII" << endl;
  vtk<< "\n";
  vtk << "DATASET UNSTRUCTURED_GRID" << endl;
  vtk << "POINTS " << nb_points << " DOUBLE" << endl;
    
  //Sortie des points
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=(P->faces).begin();F!=(P->faces).end();F++) {
	for(std::vector<int>::iterator V=faces[*F].vertex.begin();V!=faces[*F].vertex.end();V++){
	  vtk << P->mvt_t(vertex[*V].pos) << endl;
	}
      }
    }
  }
  vtk << endl;
    
  //Sortie des faces
  int point_tmp=0;
  vtk << "CELLS " << nb_faces << " " << size << endl;
  int compteur_vertex = 0;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++) {
	vtk << faces[*F].vertex.size();
	for(int k=0 ; k<faces[*F].vertex.size() ; k++)
	  vtk << " " << compteur_vertex + k; //(faces[*F].vertex)[k];
	vtk << endl;
	compteur_vertex += faces[*F].vertex.size();
      }
      //vtk << endl;
    }
  }
  vtk << endl;
  vtk << "CELL_TYPES " << nb_faces << endl;
  for(int i=0;i<nb_faces;i++){
    //vtk << 9 << endl; //Pour Quad
    vtk << 5 << endl; //Pour triangle
  }
  vtk << "\n";
  vtk << "CELL_DATA " << nb_faces << endl;
  //Deplacement au centre de la face
  vtk << "VECTORS deplacement double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	vtk << faces[*F].I_Dx << endl;
      }
    }
  }
  vtk << "\n";
  //Vitesse au centre de la face
  vtk << "VECTORS vitesse double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	vtk << faces[*F].u << endl;
      }
    }
  }
  vtk << "\n";
  //Forces au centre de la face
  vtk << "VECTORS Force double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	vtk << faces[*F].F << endl;
      }
    }
  }
  vtk << "\n";
  //Contrainte normale sur la face
  vtk << "VECTORS contrainte_normale double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	vtk << P->contrainte * faces[*F].normale << endl;
      }
    }
  }
  vtk << "\n";
  //Conditions limites sur la face
  vtk << "SCALARS BC double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	vtk << faces[*F].BC << endl;
      }
    }
  }
  vtk << "\n";
  //Indicateur de splitting sur la face
  vtk << "SCALARS split double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	vtk << faces[*F].split << endl;
      }
    }
  }
  vtk << "\n";
  //Indicateur de l'utilisation de voisins_face sur la face
  vtk << "SCALARS face_pb double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	if(faces[*F].face_pb)
	  vtk << 1. << endl;
	else
	  vtk << 0. << endl;
	  //vtk << faces[*F].face_pb << endl;
      }
    }
  }
  vtk << "\n";
  //Normale
  vtk << "VECTORS normale double" << endl;
  //vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	vtk << faces[*F].normale << endl;
      }
    }
  }
  vtk << "\n";
  //Min des coefficients de reconstruction
  /*vtk << "SCALARS c_min double 1" << endl;
  vtk << "LOOKUP_TABLE default" << endl;
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    if(not(P->split)){
      for(std::vector<int>::iterator F=P->faces.begin();F!=P->faces.end();F++){
	if((faces[*F].c_reconstruction.size()!=4 || faces[*F].reconstruction.size()!=4) && faces[*F].BC==0){
	  cout << "Face=" << *F << " Taille reconstruction=" << faces[*F].c_reconstruction.size() << " " << faces[*F].reconstruction.size() << " split=" << faces[*F].split << " BC=" << faces[*F].BC << endl;
	  getchar();
	}
	if(faces[*F].BC==0){
	  vtk << min(min(faces[*F].c_reconstruction[0],faces[*F].c_reconstruction[1]),min(faces[*F].c_reconstruction[2],faces[*F].c_reconstruction[3])) << endl;
	} else {
	  vtk << 1 << endl;
	}
      }
    }
  }
  vtk << "\n";*/
  vtk.close();
}

#endif
