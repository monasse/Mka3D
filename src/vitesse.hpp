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

#include "geometry.hpp"
#ifndef VITESSE_HPP
#define VITESSE_HPP

//Initial velocity of the solid particles
Vector_3 velocity(const Point_3 &p);


//Initial angular velocity of the solid particles
Vector_3 omega(const Point_3 &p);

Vector_3 velocity_BC(const Point_3 &p, const double& t, const double& T, const Vector_3& Dx);

Vector_3 displacement_BC(const Point_3 &p, const Vector_3 &Dx, const double& t, const double& T);

double velocity_BC_bis(const Point_3 &p, const double& t, const double& T, const Vector_3& u);


#endif
