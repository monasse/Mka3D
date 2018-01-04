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
#include "face.hpp"
#include "vertex.hpp"
#include "geometry.hpp"
#ifndef FACE_CPP
#define FACE_CPP

Face::Face()
{
  centre = Point_3(0.,0.,0.);
  normale = Vector_3(1.,0.,0.);
  voisin = -1;
  D0 = 1.;
}

Face::Face(const double& surface)
{
  centre = Point_3(0.,0.,0.);
  normale = Vector_3(1.,0.,0.);
  voisin = -1;
  D0 = 1.;
  S = surface;
}

/*Face::Face(const std::vector<Vertex> & v, const int& part)
{
  std::vector<Point_3> points;
  for(int i=0; i<v.size(); i++){
    vertex.push_back(v[i]);
    points.push_back(v[i].pos);
  }
  centre = centroid(points.begin(),points.end());
  normale = orthogonal_vector(points[0],points[1],points[2]);
  double norm = sqrt((normale.squared_length()));
  normale = normale*1./norm;
  voisin = part;
  D0 = 1.;
}*/

/*Face::Face(const std::vector<Vertex> & v, const int& part, const double& dist)
{
  std::vector<Point_3> points;
  for(int i=0; i<v.size(); i++){
    vertex.push_back(v[i]);
    points.push_back(v[i].pos);
  }
  centre = centroid(points.begin(),points.end());
  normale = orthogonal_vector(points[0],points[1],points[2]);
  double norm = sqrt((normale.squared_length()));
  normale = normale*1./norm;
  voisin = part;
  D0 = dist;
}*/

Face & Face:: operator=(const Face &F){
  assert(this != &F);
  centre = F.centre;
  normale = F.normale;
  voisin = F.voisin;
  D0  = F.D0; 
  Is = F.Is; 
  It = F.It; 
  s = F.s; 
  t = F.t;
  vertex.resize(F.vertex.size());
  for(int i= 0; i<F.vertex.size(); i++){
    vertex[i] = F.vertex[i];
  }
}

bool& operator==(const Face &F) { //Compare les faces
  if(vertex[0] != F.vertex[0] && vertex[0] != F.vertex[1] && vertex[0] != F.vertex[2])
    return false;
  else {
    if(vertex[1] != F.vertex[0] && vertex[1] != F.vertex[1] && vertex[1] != F.vertex[2])
      return false;
    else {
      if(vertex[2] != F.vertex[0] && vertex[2] != F.vertex[1] && vertex[2] != F.vertex[2])
      return false;
      else
	return true;
    }
  }
}

void Face::comp_normal(const Point_3& ext) {
  normale = orthogonal_vector(vertex[0],vertex[1],vertex[2]);
  double norm = sqrt((normale.squared_length()));
  normale = normale / norm;

  if(Vector_3(vertex[0], ext) * normale < 0.)
    normale = -1. * normale;
}

void Face::surf(const Vertex &v1, const Vertex &v2, const Vertex &v3) {
  std::vector<Vertex> aux;
  aux.push_back(v1);
  aux.push_back(v2);
  aux.push_back(v3);
  centre = centroid(aux.begin(),aux.end());
  S = 1./2.* sqrt(cross_product(Vector_3(v1.pos,v2.pos),Vector_3(v1.pos,v3.pos)).squared_length());

}

#endif
