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

#include <math.h>
#include <iostream>
#include "geometry.hpp"

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

//Boundary velocities of the solid particles
/*Vector_3 velocity_BC(const Point_3 &p, const double& t, const double& T, const Vector_3& Dx) {
  return Vector_3(0,0,-1.); // * t / T; //En m.s^-1
  }*/

  
Vector_3 velocity_BC(const Point_3 &p, const double& t, const double& T, const Vector_3& Dx)
{
  double T_p = 1.;
  double pos_x = p.x() + Dx.x();
  double pos_y = p.y() + Dx.y();
  double pos_z = p.z() + Dx.z();
  //return Vector_3(0,0,0);

  double alpha_pt = 3.1416 / 180. * 20. / T_p; //Rotation de 20° sur [0, T]
  double r = sqrt((pos_y)*(pos_y) + (pos_x)*(pos_x));
  double theta = 0.; //atan((p.y() - 0.5) / (p.x() - 0.5)); //0.;

  //Il faut écrire le vecteur e_theta avec la position actuelle (et pas initiale) de la particule !!!!
  if(pos_x <= 0. && pos_y < 0.) {
    theta = atan((-pos_y) / (-pos_x)) ;
    return r * Vector_3(sin(theta), -cos(theta), 0.) * alpha_pt;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
  }
  else if(pos_x <= 0. && pos_y > 0.) {
    theta = atan((pos_y) / (-pos_x));
    return r * Vector_3(-sin(theta), -cos(theta), 0.) * alpha_pt;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
  }
  else if(pos_x >= 0. && pos_y < 0.) {
    theta = atan((-pos_y) / (pos_x)) ;
    return r * Vector_3(sin(theta), cos(theta), 0.) * alpha_pt;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
  }
  else if(pos_x >= 0. && pos_y > 0.) {
    theta = atan((pos_y) / (pos_x));
    return r * Vector_3(-sin(theta), cos(theta), 0.) * alpha_pt;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
  }
  else
    return Vector_3(0,0,0); //Point milieu du cylindre donc bouge pas.
      
}   


//Boundary velocities of the solid particles
Vector_3 displacement_BC(const Point_3 &p, const Vector_3 &Dx, const double& t, const double& T)
{
  /*if(t < 1. * pow(10., -8.)) {
    if(p.z() <= 0.2) {
      return Vector_3(0., 0., -0.001);
    }
    else if(p.z() >= 2.8)
      return Vector_3(0., 0., 0.001);
  }
  else
  return Dx;*/
  if(p.z() <= 0.01) {
    return Vector_3(0., 0., -0.01 * t);
  }
  else if(p.z() >= 2.99)
    return Vector_3(0., 0., 0.);
  else
    return Dx;
  //Torsion
  /*double T_p = 10.;
  double pos_x = p.x() + Dx.x();
  double pos_y = p.y() + Dx.y();
  double pos_z = p.z() + Dx.z();
  //return Vector_3(0,0,0);

  double alpha_max = 0.5; //2. deg
  double alpha = 3.1416 / 180. * alpha_max / T_p * t; //Rotation de 2° sur [0, T]
  double r = sqrt((pos_y)*(pos_y) + (pos_x)*(pos_x));
  double theta = 0.; //atan((p.y() - 0.5) / (p.x() - 0.5)); //0.;

  //Deplacement imposé sur bord
  if(pos_x <= 0. && pos_y < 0. && p.z() >= 4.8) {
    theta = atan((-pos_y) / (-pos_x)) ;
    return r * Vector_3(sin(theta), -cos(theta), 0.) * alpha;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
  }
  else if(pos_x <= 0. && pos_y > 0. && p.z() >= 4.8) {
    theta = atan((pos_y) / (-pos_x));
    return r * Vector_3(-sin(theta), -cos(theta), 0.) * alpha;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
  }
  else if(pos_x >= 0. && pos_y < 0. && p.z() >= 4.8) {
    theta = atan((-pos_y) / (pos_x)) ;
    return r * Vector_3(sin(theta), cos(theta), 0.) * alpha;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
  }
  else if(pos_x >= 0. && pos_y > 0. && p.z() >= 4.8) {
    theta = atan((pos_y) / (pos_x));
    return r * Vector_3(-sin(theta), cos(theta), 0.) * alpha;// + Vector_3(0.5, 0.5, 0.); //En m.s^-1 //Origine au milieu du cylindre
  }
  else
  return Vector_3(0.,0.,0.); //Point milieu du cylindre donc bouge pas. Ou encore face en bas */
}

double displacement_BC_bis(const Point_3 &p, const Vector_3 &Dx, const double& t, const double& T)
{
  //if(t < 1. * pow(10., -8.)) {
  /*if(p.z() <= 0.2) {
    return -0.001;
  }
  else if(p.z() >= 2.8)
  return 0.001;*/
  /*else
    return Dx.vec[2];*/
  /*if(p.z() <= 0.0005) { //-
    return 0.00005 * t; //Faire par pallier comme disait Alexandre ?
  }
  else if(p.z() >= 0.0025)
  return 0.;*/
  if(p.z() <= 0.01) {
    return 0.001 * t; //-0.01
  }
  else if(p.z() >= 0.029)
    return 0.;
}

double velocity_BC_bis(const Point_3 &p, const double& t, const double& T, const Vector_3& Dx, const Vector_3& u, const int& BC) {
  //Chargement linéaire en traction
  /*if(t < 1. * pow(10., -8.)) {
    if(p.z() <= 0.2) {
      //if( t < 0.8 * T)
      return -0.05; // * t / T; //En m.s^-1
    }
    else if(p.z() >= 2.8)
      return 0.05;
    else
      return 0.;
  }
  else
  return u[2];*/
  if(p.z() <= 0.2 && BC == 1) // && p.x() <= 0.2 && p.y() <= 0.2)
    return -0.05;
  else if(p.z() >= 2.8 && BC == 1)
      return 0.;
  else
    return u[2];
}
