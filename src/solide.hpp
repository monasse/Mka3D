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
 *  \file solide.hpp
 *  \brief D&eacute;finition des classes d&eacute;crivant le solide.
 */

#ifndef SOLIDE_HPP
#define SOLIDE_HPP

#include "geometry.hpp"
#include "particule.hpp"
#include "face.hpp"
#include "vertex.hpp"
#include <vector>

class Particule;

//! Definition de la classe Solide  
class Solide
{
	
public:
  
  Solide();//:solide(std::vector<Particule>(1)){}
  Solide(const double& E, const double& nu, const double& B1, const double& n1, const double& A1, const double& H1);
  Solide(const std::vector<Particule> & Part);
  ~Solide();
  Solide & operator=(const Solide &S); // opérateur = surcharge pour l'affectation
  //void Affiche();  //fonction auxilaire utile pour les test
  int size(){
	return solide.size();
  }
  void Impression(const int &n);
  void Init(const char* s1, const char* s2, const char* s3, const bool& rep, const int& numrep, const double& rho);
  void Init(const char* s, const bool& rep, const int& numrep, const double& rho);
  void Solve_position(const double &dt, const bool &flag_2d, const double& t, const double& T);
  //void stock_def_plastique(const double &dt);
  void Solve_vitesse(const double &dt, const bool &flag_2d, const double& Amort, const double& t, const double& T);
  void Forces(const int &N_dim, const double& dt, const double& t, const double& T);
  void Forces_internes(const double& dt, const double& t, const double& T);
  void stresses(const double& t, const double& T);
  void update_triangles();
  const double Energie();
  const double Energie_potentielle();
  const double Energie_cinetique();
  double pas_temps(const double &t, const double &T, const double &cfls, const double &E, const double &nu, const double &rhos);
  bool voisins_face(int num_face); //Renvoie faux si trouve pas voisins pour reconstruction du gradient dans la face
  bool face_existe(Face f); //Renvoie vraie si la face testée est déjà das faces
  Vector_3 trouve_coord_bary(Point_3 part_1, Point_3 part_2, Point_3 voisin1, Point_3 voisin2, Point_3 centre_face);
  void reconstruction_faces_neumann(std::vector<int> num_faces, const Matrix& contrainte, const double& t, const double& V, const double& T);
  void splitting_elements(const int& num_part, const double& rho);
  void taille_maillage(); //Calcul le h du maillage comme max de tous les h des particules

  
  // private :
  std::vector<Vertex> vertex;
  std::vector<Face> faces;
  std::vector<Particule> solide; //Particules du maillage
  double h; //Taille du maillage. Max de la taille des particules

  double lambda; //Premier coeff de lamé
  double mu; //Second coefficient de lamé
  double A; //Limite élastique initiale
  double B; //Ecrouissage JC
  double n; //Ecrouissage JC
  double H; //Ecrouissage linéaire
};


#endif
