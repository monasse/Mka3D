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

#ifndef FACE_HPP
#define FACE_HPP

#include "geometry.hpp"
#include "solide.hpp"
#include "vertex.hpp"
#include <vector>

class Solide;

class Face
{
public:
  //Méthodes
  Face();//:vertex(std::vector<Vertex>(1)){}
  Face & operator=(const Face &F); //opérateur d'affectation
  int size(){
	return vertex.size();
  }
  void comp_quantities(Solide* Sol); //Computes the outward normal, the surface and the barycentre of the face
  void solve_position(const double &dt, const double& t, const double& T);
  void solve_vitesse(const double &dt, const double& t, const double& T);
  void solve_vitesse_MEMM(const double &dt, const double& t, const double& T);
  void test_fissuration();

  //Attributs
  Point_3 centre; //!< Centre de la face
  Vector_3 normale; //!< Normale sortante à la face
  Vector_3 vec_tangent_1; //Vecteur tangent 1 à la face
  Vector_3 vec_tangent_2; //Dernier vecteur tangent à la face formant une base orthonormale
  double S; //Surface de la face
  std::vector<int> vertex; //Les sommets de la face.
  int id; //Numéro de la face
  std::vector<int> voisins; //Donne le numéro des 2 voisins de la face
  std::vector<Vector_3> vitesse; //Donne les 2 vitesses pour le DDL associé à chaque particule
  //Sert pour calculer le critère de fissuration
  std::vector<int> reconstruction; //Donne le numéro des 4 particules pour la reconstruction sur la face
  std::vector<double> c_reconstruction; //Coordonnées barycentriques des centres des particules pour calcul du gradient
  double m; //Masse du DDL
  std::vector<double> masses; //Masses associées à chaque diamant dans les 2 particules qui partagent la face
  std::vector<Vector_3> Forces; //Forces sur chacun des 2 ddl...
  Vector_3 I_Dx; //Dx calculé sur la face par interpolation avec valeurs des particules du tétra
  std::vector<Vector_3> Dx; //Déplacement de cheque côté après fissuration
  Vector_3 I_Dx_prev; //Dx au pas de temps précédetn. Pour intégration MEMM
  Vector_3 u; //Vitesse des deux diamants associés
  Vector_3 u_prev; //Idem au pas de temps précédent
  Vector_3 u_prev2; //Pour intégration MEMM
  Vector_3 F; //Resultante des forces sur le DDL
  int BC; //Vaut -1 si particule au bord et peut valoir 1,2,etc... selon la condition de bord

  int type; //Pour savoir si triangle ou quad. Utilise notation gmsh
  double D0; //Distance centre face. Sera utile pour CFL...
  double h; //Diamètre de la face

  bool fissure; //Pour savoir si face a fissurée ou non

  //Booléen pour savoir si face a été splitée et doit donc être exclue des calculs
  bool split;
  //Pour indiquer s'il y a un probleme de reconstruction
  bool face_pb;
};

bool operator==(const Face &F1, const Face &F2); //Compare les faces

#endif
