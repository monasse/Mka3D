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
Vector_3 velocity_BC(const Point_3 &p, const double& t, const double& T, const Vector_3& Dx)
{
  double T_p = 0.1;
  double pos_x = p.x() + Dx.x();
  double pos_y = p.y() + Dx.y();
  double pos_z = p.z() + Dx.z();
  //return Vector_3(0,0,0);

  //Chargement linéaire en traction
  if(p.z() <= 0.2)
    return Vector_3(0,0,-0.01); // * t / T; //En m.s^-1
  else if(p.z() >= 2.8)
    return Vector_3(0,0,0);

  /*if(pos_z <= 0.1) { //Vitesse en BC...
    double alpha_pt = 3.1416 / 180. * 20. / T_p; //Rotation de 20° sur [0, T]
    double r = sqrt((pos_y-0.5)*(pos_y-0.5) + (pos_x-0.5)*(pos_x-0.5));
    double theta = 0.; //atan((p.y() - 0.5) / (p.x() - 0.5)); //0.;

    //Il faut écrire le vecteur e_theta avec la position actuelle (et pas initiale) de la particule !!!!
    if(pos_x <= 0.5 && pos_y < 0.5) {
      theta = atan((0.5 - pos_y) / (0.5 - pos_x)) ;
      return r * Vector_3(sin(theta), -cos(theta), 0.) * alpha_pt;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
    }
    else if(pos_x <= 0.5 && pos_y > 0.5) {
      theta = atan((pos_y - 0.5) / (0.5 - pos_x));
      return r * Vector_3(-sin(theta), -cos(theta), 0.) * alpha_pt;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
    }
    else if(pos_x >= 0.5 && pos_y < 0.5) {
      theta = atan((0.5 - pos_y) / (pos_x - 0.5)) ;
      return r * Vector_3(sin(theta), cos(theta), 0.) * alpha_pt;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
    }
    else if(pos_x >= 0.5 && pos_y > 0.5) {
      theta = atan((pos_y - 0.5) / (pos_x - 0.5));
      return r * Vector_3(-sin(theta), cos(theta), 0.) * alpha_pt;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
    }
    else
      return Vector_3(0,0,0); //Point milieu du cylindre donc bouge pas.
      
  }
  else if(pos_z >= 2.9)
  return Vector_3(0,0,0);*/
}

//Boundary velocities of the solid particles
Vector_3 displacement_BC(const Point_3 &p, const Vector_3 &Dx, const double& t, const double& T)
{
  double T_p = 0.001;
  
  //Chargement linéaire en traction
  /*if(p.z() <= 0.2)
    return Vector_3(-0.1,0,0) * t / T_p; //En m.
  else
  return Dx;*/

  //Chargement torsion
  /*if(p.z() <= 0.1) { //Déplacement en BC...
    double alpha_pt = 3.1416 / 360. * 45 / T; // * t / T; //Rotation de 45° sur [0, T]
    double r = sqrt((p.y()-0.5)*(p.y()-0.5) + (p.x()-0.5)*(p.x()-0.5));
    double theta = 0.;

    if(p.x() >= 0.5)
      theta = atan((p.y() - 0.5) / (p.x() - 0.5));
    else //if(p.x() >= 0.5 && p.y() >= 0.5)
      theta = 3.1416 - atan((p.y() - 0.5) / (p.x() - 0.5));
    return r * Vector_3(-sin(alpha_pt), cos(alpha_pt), 0.); //En m
  }
 else */
    return Dx;
}

double velocity_BC_bis(const Point_3 &p, const double& t, const double& T, const Vector_3& Dx, const Vector_3& u) {
  //Chargement linéaire en traction
  if(p.z() <= 0.2) {
    if( t < 0.8 * T)
      return -0.05; // * t / T; //En m.s^-1
    else
      return u[2];
  }
  else if(p.z() >= 2.8)
    return 0.;
  else
    return u[2];
}
