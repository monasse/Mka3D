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

  //Attributs
  Point_3 centre; //!< Centre de la face
  Vector_3 normale; //!< Normale sortante à la face
  double m; //Masse du diamant associé à la face
  double S; //Surface de la face
  std::vector<int> vertex; //Les sommets de la face.
  int id; //Numéro de la face
  std::vector<int> voisins; //Donne le numéro des 2 voisins de la face puis celui des 2 autres particules pour avoir le tétra associé à la face et calculer le gradient
  std::vector<int> reconstruction; //Donne le numéro des 4 particules pour la reconstruction sur la face
  std::vector<double> c_reconstruction; //Coordonnées barycentriques des centres des particules pour calcul du gradient
  Vector_3 I_Dx; //Dx calculé sur la face par interpolation avec valeurs des particules du tétra
  Vector_3 I_u; //u calculé sur la face. Interpolation dans bulk. Intégration forces au bord
  int BC; //Vaut -1 si particule au bord et peut valoir 1,2,etc... selon la condition de bord
  Vector_3 Fi; //Forces sur la face

  int type; //Pour savoir si triangle ou quad. Utilise notation gmsh
  double D0; //Distance centre face. Sera utile pour CFL...
};

bool operator==(const Face &F1, const Face &F2); //Compare les faces

#endif
