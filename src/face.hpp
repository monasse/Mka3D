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
#include <map>

class Face
{
public:
  Face();//:vertex(std::vector<Vertex>(1)){}
  Face(const double& surface);
  //Face(const std::vector<Vertex> & v, const int &part);
  //Face(const std::vector<Vertex> & v, const int &part, const double &dist);
  Face & operator=(const  Face &F); // op&eacute;rateur = surcharge pour l'affectation
  int size(){
	return vertex.size();
  }
  //void compFaceIntegrals(double &Fa, double &Fb, double &Fc, double &Faa, double &Fbb, double &Fcc, double &Faaa, double &Fbbb, double &Fccc, double &Faab, double &Fbbc, double &Fcca, const double& na, const double& nb, const double& nc, const int& a, const int& b, const int& c);
  //void compProjectionIntegrals(double &P1, double &Pa, double &Pb, double &Paa, double &Pab, double &Pbb, double &Paaa, double &Paab, double &Pabb, double &Pbbb, const int &a, const int &b, const int &c);
  //void Inertie(const Particule& part);
  Point_3 centre; //!< Centre de la face
  Vector_3 normale; //!< Normale sortante &agrave; la face
  double S; //Surface de la face
  double Is; //!< Premier moment d'inertie de la face
  double It; //!< Second moment d'inertie de la face
  Vector_3 s; //!< Vecteur selon le premier axe principal d'inertie de la face
  Vector_3 t; //!< Vecteur selon le second axe principal d'inertie de la face
  int nb_vertex;
  std::vector<int> vertex; //!< Les sommets de la face
  int voisin; //!< Le num&eacute;ro de la particule voisine. -1 si le voisin est le fluide
  double D0; //!< Distance &agrave; l'&eacute;quilibre avec la particule voisine

};

#endif
