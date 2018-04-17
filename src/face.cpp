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

#include <iostream>
#include "face.hpp"
#include "vertex.hpp"
#include "geometry.hpp"
#ifndef FACE_CPP
#define FACE_CPP

Face::Face() : I_Dx(), I_u()
{
  centre = Point_3(0.,0.,0.);
  normale = Vector_3(1.,0.,0.);
  S = 0.;
  m = 0.;
}

Face & Face:: operator=(const Face &F){
  assert(this != &F);
  centre = F.centre;
  normale = F.normale;
  vertex.resize(F.vertex.size());
  for(int i= 0; i<F.vertex.size(); i++){
    vertex[i] = F.vertex[i];
  }
  for(int i= 0; i<F.voisins.size(); i++){
    voisins[i] = F.voisins[i];
  }
  for(int i= 0; i<F.c_voisins.size(); i++){
    c_voisins[i] = F.c_voisins[i];
  }
}

bool operator==(const Face &F1, const Face &F2) { //Compare les faces
  if(F1.vertex[0] != F2.vertex[0] && F1.vertex[0] != F2.vertex[1] && F1.vertex[0] != F2.vertex[2])
    return false;
  else {
    if(F1.vertex[1] != F2.vertex[0] && F1.vertex[1] != F2.vertex[1] && F1.vertex[1] != F2.vertex[2])
      return false;
    else {
      if(F1.vertex[2] != F2.vertex[0] && F1.vertex[2] != F2.vertex[1] && F1.vertex[2] != F2.vertex[2])
      return false;
      else
	return true;
    }
  }
}

void Face::comp_quantities(const Point_3 &v1, const Point_3 &v2, const Point_3 &v3) { //, const Point_3& ext) {
  std::vector<Point_3> aux;
  aux.push_back(v1);
  aux.push_back(v2);
  aux.push_back(v3);
  centre = centroid(aux.begin(),aux.end());
  S = 1./2. * sqrt(cross_product(Vector_3(v1,v2),Vector_3(v1,v3)).squared_length());
  normale = orthogonal_vector(aux[0], aux[1], aux[2]);
  double norm = sqrt((normale.squared_length()));
  normale = normale / norm;
  D0 = 100000000000.;
  /*if(Vector_3(aux[0], ext) * normale < 0.)
    normale = -1. * normale;
  double aux2 = Vector_3(centre, v1) * normale;
  pt_face = centre + aux2 * normale;*/
}

void Face::solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  /*Dxprev = Dx;
  err_Dx = err_Dx + u * dt;
  Dx = Dx+ err_Dx;
  err_Dx = err_Dx + (Dxprev - Dx); //Version compensation de l'erreur de sommation
  */

  I_Dx = I_Dx + I_u * dt;
  /*if(I_Dx.squared_length() > 0.0000001)
    cout << "Pos face : " << I_Dx << endl;*/
}

void Face::solve_vitesse(const double& dt, const bool& flag_2d, const double& t, const double& T){
  /*u_prev = u;
  err_u= err_u + Fi * dt / m;
  u = u + err_u; //*Amort; // + velocity_BC(x0, t, T, Dx); //Conditions aux limites en vitesse ajoutées ici
  err_u = err_u + (u_prev - u); //Version compensation erreur sommation
  */
  I_u = I_u + Fi / m * dt;
}

#endif
