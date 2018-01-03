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

#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <map>
#include "geometry.hpp"




//! Definition de la classe Vertex
class Vertex 
{
public:
  Vertex();
  Vertex(const Point_3 &p);
  Vertex(const Point_3 &p, const std::vector<int> & parts);
  Vertex & operator=(const  Vertex &V); // opérateur = surcharge pour l'affectation
  Point_3 pos; //!< Coordonnées du sommet
  int num;//!< Numéro du point dans le maillage de construction
  
  int size(){
	return particules.size();
  }
  std::vector<int> particules; //!< Vecteur de particules auxquelles \a pos appartient 
};

#endif
