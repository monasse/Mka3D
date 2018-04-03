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
#include "particule.hpp"
#include "face.hpp"
#include "vertex.hpp"
#include "geometry.hpp"
#include "vitesse.hpp"
#include "forces_ext.hpp"
#include <iostream>
#include <string>
#ifndef PARTICULE_CPP
#define PARTICULE_CPP

//using namespace std;

inline double signe(const double &x)
{
  return (x < 0.) ? -1. : 1. ;
}

Particule::Particule(const int& Id):discrete_gradient(), contrainte(), epsilon_p(), vertices(), x0(), err_u(), err_Dx() {
  id = Id;
  def_plas_cumulee = 0.; //Déformation plastique cumulée du lien
  seuil_elas = 0.;
  fixe = 0;
  BC = 0;
}

Particule::Particule():discrete_gradient(), contrainte(), epsilon_p(), vertices(), x0(), err_u(), err_Dx()
{
  id = 0;
  def_plas_cumulee = 0.; //Déformation plastique cumulée du lien
  seuil_elas = 0.;
  fixe = 0;
  BC = 0;
}


Particule::~Particule(){
}

Particule & Particule:: operator=(const Particule &P){
  assert(this != &P);	
  faces = P.faces;
  fixe = P.fixe;
  m  = P.m; 
  V = P.V; 
  x0 = P.x0;
  Dx = P.Dx;
  Dxprev = P.Dxprev;
  Fi = P.Fi;
  u = P.u;
  u_prev = P.u_prev;
  mvt_t = P.mvt_t;
  mvt_tprev = P.mvt_tprev;
}

void Particule::solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  Dxprev = Dx;
  err_Dx = err_Dx + u * dt;
  Dx = Dx+ err_Dx;
  err_Dx = err_Dx + (Dxprev - Dx); //Version compensation de l'erreur de sommation
  //Dx = x0.z() * x0.z() / 9. * 4 * Vector_3(0., 0., 1.);
  //Dx = x0.z() /  3. * 4 * Vector_3(0., 0., 1.);

  //Mise a jour de la transformation donnant le mouvement de la particule
  mvt_tprev = mvt_t;
  //Aff_transformation_3 rotation(rot[0][0],rot[1][0],rot[2][0],rot[0][1],rot[1][1],rot[2][1],rot[0][2],rot[1][2],rot[2][2]);
  Aff_transformation_3 translation(Vector_3(Point_3(0.,0.,0.),x0)+Dx);
  Aff_transformation_3 translation_inv(Vector_3(x0,Point_3(0.,0.,0.)));
  mvt_t = translation*(/*rotation*/translation_inv);
	//cout<<"position du centre de la particule "<<x0+Dx<<endl;
}

void Particule::solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  u_prev = u;
  err_u = err_u + Fi * dt / m;
  u = u + err_u; //*Amort; // + velocity_BC(x0, t, T, Dx); //Conditions aux limites en vitesse ajoutées ici
  err_u = err_u + (u_prev - u); //Version compensation erreur sommation
  u.vec[2] = velocity_BC_bis(x0, t, T, Dx, u);
}

void Particule::barycentre() {
  x0 = centroid(vertices.begin(), vertices.end());
}

void Particule::volume() {
  V = abs(cross_product(Vector_3(vertices[0],vertices[1]),Vector_3(vertices[0],vertices[2]))*Vector_3(vertices[0],vertices[3]))/6.;
}

#endif
