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
#include "vertex.hpp"
#include "geometry.hpp"
#ifndef VERTEX_CPP
#define VERTEX_CPP

Vertex::Vertex()
{
  pos = Point_3(0.,0.,0.);
  num = 0;
}

Vertex::Vertex(const Point_3 &p) {
  pos = p;
}

Vertex::Vertex(const Point_3& p, const std::vector<int> & parts)
{
  pos = p;
  for(int i=0; i<parts.size(); i++){
    particules.push_back(parts[i]);
  }
}

Vertex & Vertex:: operator=(const Vertex &V){
	
  assert(this != &V);
  pos = V.pos;
  num = V.num;
  particules.resize(V.particules.size());
  for(int i=0; i<V.particules.size(); i++){
    particules[i]= V.particules[i];
  }
}

#endif
