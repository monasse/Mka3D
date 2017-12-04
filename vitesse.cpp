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
#include <math.h>

//Initial velocity of the solid particles
Vector_3 velocity(const Point_3 &p)
{
  return Vector_3(0,0,0);
  //return Vector_3(0.01*(p.x()-1.),0,0);
  /*if(p.x() <= 1.)
     return Vector_3(-20.,0,0);
  else if(p.x() >= 19.)
  return Vector_3(0,0,0); */
}

//Initial angular velocity of the solid particles
Vector_3 omega(const Point_3 &p)
{
  return Vector_3(0,0,0);
  //return Vector_3(0.1,0,0);
  
}

//Boundary velocities of the solid particles
Vector_3 velocity_BC(const Point_3 &p, const double& t, const double& T)
{
  return Vector_3(0,0,0);

  //Chargement linéaire en traction
  /*if(p.x() <= 1.)
    return Vector_3(-5.,0,0); // * t / T; //En m.s^-1
  else if(p.x() >= 4.)
  return Vector_3(0,0,0);*/

  /*if(p.x() <= 1.) { //Vitesse en BC...
    double alpha_pt = 3.1416 / 18. / T; //Rotation de 20° sur [0, T]
    //double r = sqrt((p.y()-0.5)*(p.y()-0.5) + (p.z()-0.5)*(p.z()-0.5));
    double theta = atan((p.z() - 0.5) / (p.y() - 0.5)) ;
    return Vector_3(0.,-sin(theta), cos(theta)) * alpha_pt * t; //En m.s^-1
  }
  else if(p.x() >= 4.)
  return Vector_3(0,0,0);*/
}

//Boundary velocities of the solid particles
Vector_3 displacement_BC(const Point_3 &p, const Vector_3 &Dx, const double& t, const double& T)
{
  //Chargement linéaire en traction
  /*if(p.x() <= 1.)
    return Vector_3(-0.1,0,0) * t / T; //En m.
  else
  return Dx; */

  //Chargement torsion
  if(p.z() >= 2.8) { //Vitesse en BC...
    double alpha_pt = 3.1416 / 360. * 45 * t / T; //Rotation de 20° sur [0, T]
    double r = sqrt((p.y()-0.5)*(p.y()-0.5) + (p.x()-0.5)*(p.x()-0.5));

    if(p.x() >= 0.5)
      theta = atan((p.y() - 0.5) / (p.x() - 0.5));
    else //if(p.x() >= 0.5 && p.y() >= 0.5)
      theta = 3.1416 - atan((p.y() - 0.5) / (p.x() - 0.5));
	//else if(p.x() >= 0.5 && p.y() <= 0.5)
    return r * Vector_3(cos(theta), sin(theta), 0.) * alpha_pt * t; //En m.s^-1
  }
  else
    return Dx;
}
