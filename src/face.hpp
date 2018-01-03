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
#include <vector>

class Face
{
public:
  //Méthodes
  Face();//:vertex(std::vector<Vertex>(1)){}
  Face(const double& surface);
  //Face(const std::vector<Vertex> & v, const int &part);
  //Face(const std::vector<Vertex> & v, const int &part, const double &dist);
  Face & operator=(const  Face &F); //opérateur d'affectation
  int size(){
	return vertex.size();
  }
  void comp_normal(const Point_3& ext); //Calcul une normale à la face
  bool& operator==(const Face &F); //Compare les faces

  //Attributs
  //Point_3 centre; //!< Centre de la face
  Vector_3 normale; //!< Normale sortante &agrave; la face
  double S; //Surface de la face
  double Is; //!< Premier moment d'inertie de la face
  double It; //!< Second moment d'inertie de la face
  Vector_3 s; //!< Vecteur selon le premier axe principal d'inertie de la face
  Vector_3 t; //!< Vecteur selon le second axe principal d'inertie de la face
  int nb_vertex;
  std::vector<int> vertex; //!< Les sommets de la face
  //int voisin; //!< Le numéro de la particule voisine. -1 si le voisin est le fluide
  double D0; //!< Distance à l'équilibre avec la particule voisine
  std::vector<Particule *> parts; //contient les deux particules dans le lien. La normale est dans le sens de la première vers la seconde

};

#endif
