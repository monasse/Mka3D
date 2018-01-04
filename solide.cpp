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
#include "solide.hpp"
#include "geometry.hpp"
#include "forces_ext.cpp"
#include "vitesse.cpp"
#include <iostream>
#ifndef SOLIDE_CPP
#define SOLIDE_CPP



inline double signe(const double &x)
{
  return (x < 0.) ? -1. : 1. ;
}

//Function absolute value
double abs(const double &x)
{
  return max(x,-x);
}


//const double eps_relat = numeric_limits<double>::epsilon();
const double eps_relat =0.000001;
/*!
* \fn Vertex::Vertex()
* \brief Default constructor 
*/
Vertex::Vertex()
{
  pos = Point_3(0.,0.,0.);
  num = 0;
}

/*!
 *\fn Vertex::Vertex(const Point_3 p, std::vector<int> & parts)
 *\brief Constructor overload.
 *\param p: coordinates of the vertex
 *\param parts: vector of particles sharing the vertex
 */
Vertex::Vertex(const Point_3& p, const std::vector<int> & parts)
{
  pos = p;
  for(int i=0; i<parts.size(); i++){
    particules.push_back(parts[i]);
  }
}
/*!
* \fn Vertex & Vertex:: operator=(const Vertex &V)
* \brief operator =: affectation overload.
* \param V Vertex
* \return Vertex
*/
Vertex & Vertex:: operator=(const Vertex &V){
	
	assert(this != &V);
	pos = V.pos;
	num = V.num;
	particules.resize(V.particules.size());
	for(int i=0; i<V.particules.size(); i++){
		particules[i]= V.particules[i];
	}
}
/*!
* \fn Face::Face()
* \brief Default constructor. 
 */
Face::Face()
{
  centre = Point_3(0.,0.,0.);
  normale = Vector_3(1.,0.,0.);
  voisin = -1;
  D0 = 1.;
}

/*!
 *\fn Face::Face(std::vector<Vertex> & v, int part)
 *\brief Surcharge du constructeur
 *\param v vecteur de sommets 
 *\param part num&eacute;ro de la particule voisine. -1 si le voisin est le fluide 
 */
Face::Face(const std::vector<Vertex> & v, const int& part)
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
}
/*!
* \fn Face::Face(std::vector<Vertex> & v, int part, double dist)
* \brief Surcharge du constructeur.
* \param v vecteur de sommets
* \param part num&eacute;ro de la particule voisine. -1 si le voisin est le fluide 
* \param dist distance &agrave; l'&eacute;quilibre avec la particule voisine
 */
Face::Face(const std::vector<Vertex> & v, const int& part, const double& dist)
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
}
/*!
* \fn Face & Face:: operator=(const Face &F)
* \brief op&eacute;rateur =
* \param F Face
* \return Face
 */
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
/*!
* \fn void Face::Inertie()
* \brief Calcul d'inertie de la face. 
* \details 
 * Premier moment d'inertie de la face \n
 * Second moment d'inertie de la face \n
 * Vecteur selon le premier axe principal d'inertie de la face \n
 * Vecteur selon le second axe principal d'inertie de la face 
 * \return void
 */
void Face::Inertie(){
  double eps = 1e-14;//std::numeric_limits<double>::epsilon();
  //Choix initial d'un repere orthonorme de la face
  if(normale[0]!=0. || normale[1]!=0.){
    s = Vector_3(-normale[1],normale[0],0.);
    s = s/(sqrt((s.squared_length())));
    t = cross_product(normale,s);
  } else {
    s = Vector_3(0.,0.,1.);
    s = s/(sqrt((s.squared_length())));
    t = cross_product(normale,s);
  }
  //Calcul du centre de la face
  double T1 = 0.;
  double Ts = 0.;
  double Tt = 0.;
  for(int i=0;i<size();i++){
    int ip = (i+1)%(size());
    Vector_3 v1(centre,vertex[i].pos);
    Vector_3 v2(centre,vertex[ip].pos);
    T1 += 1./2.*(cross_product(v1,v2)*normale);
    Ts += 1./6.*(cross_product(v1,v2)*normale)*((v1+v2)*s);
    Tt += 1./6.*(cross_product(v1,v2)*normale)*((v1+v2)*t);
  }
  centre = centre + (Ts/T1)*s + (Tt/T1)*t;
  S = T1;
  //Calcul de la matrice d'inertie de la face dans les deux axes a l'origine centre
  double Tss = 0.;
  double Ttt = 0.;
  double Tst = 0.;
  for(int i=0;i<size();i++){
    int ip = (i+1)%(size());
    Vector_3 v1(centre,vertex[i].pos);
    Vector_3 v2(centre,vertex[ip].pos);
    double As = (v1*s);
    double At = (v1*t);
    double Bs = (v2*s);
    double Bt = (v2*t);
    Tss += 1./12.*(As*As+As*Bs+Bs*Bs);
    Ttt += 1./12.*(At*At+At*Bt+Bt*Bt);
    Tst += 1./24.*(2.*As*At+As*Bt+At*Bs+2.*Bs*Bt);
  }
  //Calcul des moments d'inertie
  double Delta = pow(Tss-Ttt,2)+4.*Tst*Tst;
  Is = (Tss+Ttt+sqrt(Delta))/2.;
  It = (Tss+Ttt-sqrt(Delta))/2.;
  //Diagonalisation
  if(abs(Tss-Ttt)>eps){
    if(abs(Tss-Is)>eps){
      Vector_3 stemp = -Tst*s+(Tss-Is)*t;
      s = stemp/(sqrt((stemp.squared_length())));
      t = cross_product(normale,s);
    } else {
      Vector_3 stemp = -Tst*t+(Ttt-Is)*s;
      s = stemp/(sqrt((stemp.squared_length())));
      t = cross_product(normale,s);
    }
  } else {
    if(abs(Tst)>eps){
      Vector_3 stemp = s+t;
      Vector_3 ttemp = -s+t;
      s = stemp/(sqrt((stemp.squared_length())));
      t = stemp/(sqrt((ttemp.squared_length())));
    }
  }
}


/*!
 * \fn Particule::Particule()
 * \brief Constructeur par d&eacute;faut. 
 */

Particule::Particule():discrete_gradient(), contrainte(), epsilon_p(), n_elas_prev()
{
  //discrete_gradient = 0.; //Gradient reconstruit par particule
  //contrainte = 0.; //Contrainte par particule
  def_plas_cumulee = 0.; //Déformation plastique cumulée du lien
  seuil_elas = 0.;
  //epsilon_p = 0.;
  
  /*min_x = 0.; 
  min_y = 0.;
  min_z = 0.;
  max_x = 1. ;
  max_y = 1.;
  max_z = 1.;*/
  bbox = Bbox(0.,0.,0.,1.,1.,1.);
  
  x0 = Point_3(0.5,0.5,0.5);
	
  cube = true;
	
  const Point_3 s1(0.,0.,0.);
  const Point_3 r1(1.,0.,0.);
  const Point_3 t1(1.,1.,0.);
  const Point_3 v1(0.,1.,0.);
	
	
  const Point_3 s2(0.,0.,1.);
  const Point_3 r2(1.,0.,1.);
  const Point_3 t2(1.,1.,1.);
  const Point_3 v2(0.,1.,1.);

  //Face 1
  std::vector<Vertex> vert1(4);
  vert1[0].pos = s1;
  vert1[0].particules.push_back(-1);
  vert1[1].pos = v1;
  vert1[1].particules.push_back(-1);
  vert1[2].pos = t1;
  vert1[2].particules.push_back(-1);
  vert1[3].pos = r1;
  vert1[3].particules.push_back(-1);
  vector<Point_3> points1(4);
  for(int i=0;i<4;i++){
    points1[i] = vert1[i].pos;
  }
  Face face1 = Face(vert1,-1,1.);
  face1.centre = centroid(points1.begin(),points1.end());
  face1.normale = orthogonal_vector(points1[0],points1[1],points1[2]);
  double norm1 = sqrt((face1.normale.squared_length()));
  face1.normale = face1.normale*1./norm1;

  //Face 2
  std::vector<Vertex> vert2(4);
  vert2[0].pos = s1;
  vert2[0].particules.push_back(-1);
  vert2[1].pos = r1;
  vert2[1].particules.push_back(-1);
  vert2[2].pos = r2;
  vert2[2].particules.push_back(-1);
  vert2[3].pos = s2;
  vert2[3].particules.push_back(-1);
  vector<Point_3> points2(4);
  for(int i=0;i<4;i++){
    points2[i] = vert2[i].pos;
  }
  Face face2 = Face(vert2,-1,1.);
  face2.centre = centroid(points2.begin(),points2.end());
  face2.normale = orthogonal_vector(points2[0],points2[1],points2[2]);
  double norm2 = sqrt((face2.normale.squared_length()));
  face2.normale = face2.normale*1./norm2;

  //Face 3
  std::vector<Vertex> vert3(4);
  vert3[0].pos = r1;
  vert3[0].particules.push_back(-1);
  vert3[1].pos = t1;
  vert3[1].particules.push_back(-1);
  vert3[2].pos = t2;
  vert3[2].particules.push_back(-1);
  vert3[3].pos = r2;
  vert3[3].particules.push_back(-1);
  vector<Point_3> points3(4);
  for(int i=0;i<4;i++){
    points3[i] = vert3[i].pos;
  }
  Face face3 = Face(vert3,-1,1.);
  face3.centre = centroid(points3.begin(),points3.end());
  face3.normale = orthogonal_vector(points3[0],points3[1],points3[2]);
  double norm3 = sqrt((face3.normale.squared_length()));
  face3.normale = face3.normale*1./norm3;

  //Face 4
  std::vector<Vertex> vert4(4);
  vert4[0].pos = t1;
  vert4[0].particules.push_back(-1);
  vert4[1].pos = v1;
  vert4[1].particules.push_back(-1);
  vert4[2].pos = v2;
  vert4[2].particules.push_back(-1);
  vert4[3].pos = t2;
  vert4[3].particules.push_back(-1);
  vector<Point_3> points4(4);
  for(int i=0;i<4;i++){
    points4[i] = vert4[i].pos;
  }
  Face face4 = Face(vert4,-1,1.);
  face4.centre = centroid(points4.begin(),points4.end());
  face4.normale = orthogonal_vector(points4[0],points4[1],points4[2]);
  double norm4 = sqrt((face4.normale.squared_length()));
  face4.normale = face4.normale*1./norm4;
	
  //Face 5
  std::vector<Vertex> vert5(4);
  vert5[0].pos = s1;
  vert5[0].particules.push_back(-1);
  vert5[1].pos = r1;
  vert5[1].particules.push_back(-1);
  vert5[2].pos = r2;
  vert5[2].particules.push_back(-1);
  vert5[3].pos = s2;
  vert5[3].particules.push_back(-1);
  vector<Point_3> points5(4);
  for(int i=0;i<4;i++){
    points5[i] = vert5[i].pos;
  }
  Face face5 = Face(vert5,-1,1.);
  face5.centre = centroid(points5.begin(),points5.end());
  face5.normale = orthogonal_vector(points5[0],points5[1],points5[2]);
  double norm5 = sqrt((face5.normale.squared_length()));
  face5.normale = face5.normale*1./norm5;
  
  //Face 6
  std::vector<Vertex> vert6(4);
  vert6[0].pos = s1;
  vert6[0].particules.push_back(-1);
  vert6[1].pos = r1;
  vert6[1].particules.push_back(-1);
  vert6[2].pos = r2;
  vert6[2].particules.push_back(-1);
  vert6[3].pos = s2;
  vert6[3].particules.push_back(-1);
  vector<Point_3> points6(4);
  for(int i=0;i<4;i++){
    points6[i] = vert6[i].pos;
  }
  Face face6 = Face(vert6,-1,1.);
  face6.centre = centroid(points6.begin(),points6.end());
  face6.normale = orthogonal_vector(points6[0],points6[1],points6[2]);
  double norm6 = sqrt((face6.normale.squared_length()));
  face6.normale = face6.normale*1./norm6;

  std::vector<Face> f(6);
  f[0] = face1;
  f[1] = face2;
  f[2] = face3;
  f[3] = face4;
  f[4] = face5;
  f[5] = face6;

  faces = f;
	
  //triangles par face	
  Triangle_3 Tri1(s1,r1,v1);
  Triangle_3 Tri2(t1,r1,v1);	
  triangles.push_back(Tri1);
  triangles.push_back(Tri2);
  normales.push_back(face1.normale);
  normales.push_back(face1.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face2
  Triangle_3 Tri5(s2,r2,v2);
  Triangle_3 Tri6(t2,r2,v2);
  triangles.push_back(Tri5);
  triangles.push_back(Tri6);
  normales.push_back(face2.normale);
  normales.push_back(face2.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face3
  Triangle_3 Tri9(s2,s1,v2);
  Triangle_3 Tri10(v1,s1,v2);
  triangles.push_back(Tri9);
  triangles.push_back(Tri10);
  normales.push_back(face3.normale);
  normales.push_back(face3.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face4	
  Triangle_3 Tri13(r2,r1,t2);
  Triangle_3 Tri14(t1,r1,t2);
  triangles.push_back(Tri13);
  triangles.push_back(Tri14);
  normales.push_back(face4.normale);
  normales.push_back(face4.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face5	
  Triangle_3 Tri17(v2,v1,t2);
  Triangle_3 Tri18(t1,v1,t2);
  triangles.push_back(Tri17);
  triangles.push_back(Tri18);
  normales.push_back(face5.normale);
  normales.push_back(face5.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face6
  Triangle_3 Tri21(s2,s1,r2);
  Triangle_3 Tri22(r1,s1,r2);
  triangles.push_back(Tri21);
  triangles.push_back(Tri22); 
  normales.push_back(face6.normale);
  normales.push_back(face6.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  Points_interface.resize(triangles.size(), std::vector<Point_3>(0));
  Triangles_interface.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Points_interface_prev.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface_prev.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface_prev.resize(triangles.size(), std::vector<std::vector<int> >(0));
	Ff = Vector_3(0.,0.,0.); Ffprev = Vector_3(0.,0.,0.); 
	Mf = Vector_3(0.,0.,0.); Mfprev = Vector_3(0.,0.,0.);
}

/**
*\fn Particule::Particule(const double x_min, const double y_min, const double z_min,  const double x_max, const double y_max,const double z_max)
*\brief Surcharge du constructeur.
* \param (x_min, y_min, z_min) : coordonn&eacute;es du sommet le plus &agrave; gauche de la particule
* \param (x_max, y_max, z_max) : coordonn&eacute;es du sommet le plus &agrave; droite de la particule
*/

Particule::Particule(const double &x_min, const double &y_min, const double &z_min, 
		     const double &x_max, const double &y_max, const double &z_max) : discrete_gradient(), contrainte(), epsilon_p(), n_elas_prev()
{
  //discrete_gradient = 0.; //Gradient reconstruit par particule
  //contrainte = 0.; //Contrainte par particule
  def_plas_cumulee = 0.; //Déformation plastique cumulée du lien
  seuil_elas = 0.;
  //epsilon_p = 0.;
  
  /*min_x = x_min; 
  min_y = y_min;
  min_z = z_min;
  max_x = x_max ;
  max_y = y_max;
  max_z = z_max;*/
  bbox = Bbox(x_min,y_min,z_min,x_max,y_max,z_max);
  
  x0 = Point_3((x_min+x_max)/2.,(y_min+y_max)/2.,(z_min+z_max)/2.);
	
  cube = true;
	
  const Point_3 s1(x_min, y_min, z_min);
  const Point_3 r1(x_max, y_min, z_min);
  const Point_3 t1(x_max, y_max, z_min);
  const Point_3 v1(x_min, y_max, z_min);
	
	
  const Point_3 s2(x_min, y_min, z_max);
  const Point_3 r2(x_max, y_min, z_max);
  const Point_3 t2(x_max, y_max, z_max);
  const Point_3 v2(x_min, y_max, z_max);
	
  //Face 1
  std::vector<Vertex> vert1(4);
  vert1[0].pos = s1;
  vert1[0].particules.push_back(-1);
  vert1[1].pos = v1;
  vert1[1].particules.push_back(-1);
  vert1[2].pos = t1;
  vert1[2].particules.push_back(-1);
  vert1[3].pos = r1;
  vert1[3].particules.push_back(-1);
  vector<Point_3> points1(4);
  for(int i=0;i<4;i++){
    points1[i] = vert1[i].pos;
  }
  Face face1 = Face(vert1,-1,1.);
  face1.centre = centroid(points1.begin(),points1.end());
  face1.normale = orthogonal_vector(points1[0],points1[1],points1[2]);
  double norm1 = sqrt((face1.normale.squared_length()));
  face1.normale = face1.normale*1./norm1;
	
  //Face 2
  std::vector<Vertex> vert2(4);
  vert2[0].pos = s1;
  vert2[0].particules.push_back(-1);
  vert2[1].pos = r1;
  vert2[1].particules.push_back(-1);
  vert2[2].pos = r2;
  vert2[2].particules.push_back(-1);
  vert2[3].pos = s2;
  vert2[3].particules.push_back(-1);
  vector<Point_3> points2(4);
  for(int i=0;i<4;i++){
    points2[i] = vert2[i].pos;
  }
  Face face2 = Face(vert2,-1,1.);
  face2.centre = centroid(points2.begin(),points2.end());
  face2.normale = orthogonal_vector(points2[0],points2[1],points2[2]);
  double norm2 = sqrt((face2.normale.squared_length()));
  face2.normale = face2.normale*1./norm2;
	
  //Face 3
  std::vector<Vertex> vert3(4);
  vert3[0].pos = r1;
  vert3[0].particules.push_back(-1);
  vert3[1].pos = t1;
  vert3[1].particules.push_back(-1);
  vert3[2].pos = t2;
  vert3[2].particules.push_back(-1);
  vert3[3].pos = r2;
  vert3[3].particules.push_back(-1);
  vector<Point_3> points3(4);
  for(int i=0;i<4;i++){
    points3[i] = vert3[i].pos;
  }
  Face face3 = Face(vert3,-1,1.);
  face3.centre = centroid(points3.begin(),points3.end());
  face3.normale = orthogonal_vector(points3[0],points3[1],points3[2]);
  double norm3 = sqrt((face3.normale.squared_length()));
  face3.normale = face3.normale*1./norm3;

  //Face 4
  std::vector<Vertex> vert4(4);
  vert4[0].pos = t1;
  vert4[0].particules.push_back(-1);
  vert4[1].pos = v1;
  vert4[1].particules.push_back(-1);
  vert4[2].pos = v2;
  vert4[2].particules.push_back(-1);
  vert4[3].pos = t2;
  vert4[3].particules.push_back(-1);
  vector<Point_3> points4(4);
  for(int i=0;i<4;i++){
    points4[i] = vert4[i].pos;
  }
  Face face4 = Face(vert4,-1,1.);
  face4.centre = centroid(points4.begin(),points4.end());
  face4.normale = orthogonal_vector(points4[0],points4[1],points4[2]);
  double norm4 = sqrt((face4.normale.squared_length()));
  face4.normale = face4.normale*1./norm4;

  //Face 5
  std::vector<Vertex> vert5(4);
  vert5[0].pos = s1;
  vert5[0].particules.push_back(-1);
  vert5[1].pos = r1;
  vert5[1].particules.push_back(-1);
  vert5[2].pos = r2;
  vert5[2].particules.push_back(-1);
  vert5[3].pos = s2;
  vert5[3].particules.push_back(-1);
  vector<Point_3> points5(4);
  for(int i=0;i<4;i++){
    points5[i] = vert5[i].pos;
  }
  Face face5 = Face(vert5,-1,1.);
  face5.centre = centroid(points5.begin(),points5.end());
  face5.normale = orthogonal_vector(points5[0],points5[1],points5[2]);
  double norm5 = sqrt((face5.normale.squared_length()));
  face5.normale = face5.normale*1./norm5;

  //Face 6
  std::vector<Vertex> vert6(4);
  vert6[0].pos = s1;
  vert6[0].particules.push_back(-1);
  vert6[1].pos = r1;
  vert6[1].particules.push_back(-1);
  vert6[2].pos = r2;
  vert6[2].particules.push_back(-1);
  vert6[3].pos = s2;
  vert6[3].particules.push_back(-1);
  vector<Point_3> points6(4);
  for(int i=0;i<4;i++){
    points6[i] = vert6[i].pos;
  }
  Face face6 = Face(vert6,-1,1.);
  face6.centre = centroid(points6.begin(),points6.end());
  face6.normale = orthogonal_vector(points6[0],points6[1],points6[2]);
  double norm6 = sqrt((face6.normale.squared_length()));
  face6.normale = face6.normale*1./norm6;

  std::vector<Face> f(6);
  f[0] = face1;
  f[1] = face2;
  f[2] = face3;
  f[3] = face4;
  f[4] = face5;
  f[5] = face6;

  faces = f;
	
  //triangles par face	
  Triangle_3 Tri1(s1,r1,v1);
  Triangle_3 Tri2(t1,r1,v1);	
  triangles.push_back(Tri1);
  triangles.push_back(Tri2);
  normales.push_back(face1.normale);
  normales.push_back(face1.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face2
  Triangle_3 Tri5(s2,r2,v2);
  Triangle_3 Tri6(t2,r2,v2);
  triangles.push_back(Tri5);
  triangles.push_back(Tri6);
  normales.push_back(face2.normale);
  normales.push_back(face2.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face3
  Triangle_3 Tri9(s2,s1,v2);
  Triangle_3 Tri10(v1,s1,v2);
  triangles.push_back(Tri9);
  triangles.push_back(Tri10);
  normales.push_back(face3.normale);
  normales.push_back(face3.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face4	
  Triangle_3 Tri13(r2,r1,t2);
  Triangle_3 Tri14(t1,r1,t2);
  triangles.push_back(Tri13);
  triangles.push_back(Tri14);
  normales.push_back(face4.normale);
  normales.push_back(face4.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face5	
  Triangle_3 Tri17(v2,v1,t2);
  Triangle_3 Tri18(t1,v1,t2);
  triangles.push_back(Tri17);
  triangles.push_back(Tri18);
  normales.push_back(face5.normale);
  normales.push_back(face5.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  //face6
  Triangle_3 Tri21(s2,s1,r2);
  Triangle_3 Tri22(r1,s1,r2);
  triangles.push_back(Tri21);
  triangles.push_back(Tri22); 
  normales.push_back(face6.normale);
  normales.push_back(face6.normale);
  fluide.push_back(true);
  fluide.push_back(true);
  Points_interface.resize(triangles.size(), std::vector<Point_3>(0));
  Triangles_interface.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Points_interface_prev.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface_prev.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface_prev.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Ff = Vector_3(0.,0.,0.); Ffprev = Vector_3(0.,0.,0.); 
	Mf = Vector_3(0.,0.,0.); Mfprev = Vector_3(0.,0.,0.);
}

/*!
* \fn Particule::Particule(Point_3 c, const double x_min, const double y_min, const double z_min, const double x_max, const double y_max,const double z_max, std::vector<Face> & F)
* \brief Surcharge du constructeur.
* \param (x_min, y_min, z_min) : coordonn&eacute;es du sommet le plus &agrave; gauche de la particule
* \param (x_max, y_max, z_max) : coordonn&eacute;es du sommet le plus &agrave; droite de la particule
* \param c Point
* \param F : Face de la Particule
 */
Particule::Particule(const Point_3& c, const double &x_min, const double& y_min, const double& z_min, 
		     const double& x_max, const double& y_max,const double& z_max, 
		     const std::vector<Face> & F) : discrete_gradient(), contrainte(), epsilon_p(), n_elas_prev()
{
  //discrete_gradient = 0.; //Gradient reconstruit par particule
  //contrainte = 0.; //Contrainte par particule
  def_plas_cumulee = 0.; //Déformation plastique cumulée du lien
  seuil_elas = 0.;
  //epsilon_p = 0.;
  
  /*min_x = x_min; 
  min_y = y_min;
  min_z = z_min;
  max_x = x_max ;
  max_y = y_max;
  max_z = z_max;*/
  bbox = Bbox(x_min,y_min,z_min,x_max,y_max,z_max);
  
  x0 = c;
	
  cube = false;

  faces = F;

  for(int i=0;i<faces.size();i++){
		
 		if(faces[i].size() == 3){
			Point_3 s,r,v;
			s = faces[i].vertex[0].pos;
			r = faces[i].vertex[1].pos;
			v = faces[i].vertex[2].pos;
			
			Vector_3 vect0(s,r);
			Vector_3 vect1(s,v);
			Vector_3 normale = cross_product(vect0,vect1);
			normale = normale*(1./sqrt((normale.squared_length())));
			if (normale*faces[i].normale > 0.){
			Triangle_3 Tri(s,r,v);
			triangles.push_back(Tri);
			normales.push_back(faces[i].normale);
			if(faces[i].voisin < 0){
				fluide.push_back(true);
			} else {
				fluide.push_back(false);
			}
			}
			else{
				Triangle_3 Tri(s,v,r);
				triangles.push_back(Tri);
				normales.push_back(faces[i].normale);
				if(faces[i].voisin < 0){
					fluide.push_back(true);
				} else {
					fluide.push_back(false);
				}
			}
			if(faces[i].voisin == -2){
				vide.push_back(true);
			} else {
				vide.push_back(false);
			}
		}
		
// 	else if(flag_2d){
// 				Point_3 s,r,v,t;
// 				
// 
// 					s = faces[i].vertex[0].pos;
// 					r = faces[i].vertex[1].pos;
// 					v = faces[i].vertex[2].pos;
// 					t = faces[i].vertex[3].pos;
// 					
// 					Triangle_3 Tri1(s,r,v);
// 					triangles.push_back(Tri1);
// 					Triangle_3 Tri2(v,t,s);
// 					triangles.push_back(Tri2);
// 					
// 					normales.push_back(faces[i].normale);
// 					normales.push_back(faces[i].normale);
// 					
// 					if(faces[i].voisin < 0){
// 			           fluide.push_back(true);
// 								 fluide.push_back(true);
// 					} else {
// 			           fluide.push_back(false);
// 								 fluide.push_back(false);
// 					}
// 					if(faces[i].voisin == -2){
// 						vide.push_back(true);
// 						vide.push_back(true);
// 					} else {
// 						vide.push_back(false);
// 						vide.push_back(false);
// 					}
// 
// 		}
		else{
			Point_3 s,r,v;
			s = faces[i].centre;
			for(int k=0;k<faces[i].size();k++){
				int kp = (k+1)%(faces[i].size());
				r = faces[i].vertex[k].pos;
				v = faces[i].vertex[kp].pos;
				Vector_3 vect0(s,r);
				Vector_3 vect1(s,v);
				
				Triangle_3 Tri(s,r,v);
				triangles.push_back(Tri);
				normales.push_back(faces[i].normale);
				
				if(faces[i].voisin < 0){
					fluide.push_back(true);
				} else {
					fluide.push_back(false);
				}
				if(faces[i].voisin == -2){
					vide.push_back(true);
				} else {
					vide.push_back(false);
				}
			}
		}
		
  }// end boucle sur les faces
  
  Points_interface.resize(triangles.size(), std::vector<Point_3>(0));
  Triangles_interface.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Points_interface_prev.resize(triangles.size(), std::vector<Point_3>(0));
	Triangles_interface_prev.resize(triangles.size(), std::vector<Triangle_3>(0));
	Position_Triangles_interface_prev.resize(triangles.size(), std::vector< std::vector<int> >(0));
	Ff = Vector_3(0.,0.,0.); Ffprev = Vector_3(0.,0.,0.); 
	Mf = Vector_3(0.,0.,0.); Mfprev = Vector_3(0.,0.,0.);
}
/*!
* \fn Particule::~Particule()
* \brief Destructeur.
 */
Particule::~Particule(){
}

/*!
* \fn Particule & Particule:: operator=(const Particule &P)
* \brief op&eacute;rateur =
* \param P Particule
* \return Particule
 */
Particule & Particule:: operator=(const Particule &P){
	
	assert(this != &P);
	
	/*min_x = P.min_x;
	min_y = P.min_y;
	min_z = P.min_z;
	max_x = P.max_x;
	max_y = P.max_y;
	max_z = P.max_z;*/
	bbox = P.bbox;
	cube  = P.cube;
	
	faces = P.faces;
	fixe = P.fixe;
	m  = P.m; 
	V = P.V; 
	Vl = P.Vl; 
	epsilon = P.epsilon; 
	for(int i=0; i<3;i++){
		I[i] = P.I[i];
		for(int j=0; j<3;j++){
			rotref[i][j] = P.rotref[i][j];
		}
	}
	
	x0 = P.x0;
	Dx = P.Dx;
	Dxprev = P.Dxprev;
	Fi = P.Fi;
	Ff = P.Ff;
	Ffprev = P.Ffprev;
	Mi = P.Mi;
	Mf = P.Mf;
	Mfprev = P.Mfprev;
	u = P.u;
	u_half = P.u_half;
	omega = P.omega;
	omega_half = P.omega_half;
	e = P.e;
	eprev = P.eprev;
	mvt_t = P.mvt_t;
	mvt_tprev = P.mvt_tprev;
	
	triangles.resize(P.triangles.size()); 
	for(int i = 0; i< P.triangles.size(); i++){
		triangles[i] = P.triangles[i];
	}
	vide.resize(P.vide.size()); 
	for(int i = 0; i< P.vide.size(); i++){
		vide[i] = P.vide[i];
	}
	triangles_prev.resize(P.triangles_prev.size()); 
	for(int i = 0; i< P.triangles_prev.size(); i++){
		triangles_prev[i] = P.triangles_prev[i];
	}
	
	normales.resize(P.normales.size()); 
	for(int i = 0; i< P.normales.size(); i++){
		normales[i] = P.normales[i];
	}
	normales_prev.resize(P.normales_prev.size()); 
	for(int i = 0; i< P.normales_prev.size(); i++){
		normales_prev[i] = P.normales_prev[i];
	}
	fluide.resize(P.fluide.size()); 
	for(int i = 0; i< P.fluide.size(); i++){
		fluide[i] = P.fluide[i];
	}
	fluide_prev.resize(P.fluide_prev.size()); 
	for(int i = 0; i< P.fluide_prev.size(); i++){
		fluide_prev[i] = P.fluide_prev[i];
	}
	
	Points_interface.resize(P.Points_interface.size(), std::vector<Point_3>(0));
	for(int i= 0; i< P.triangles.size(); i++){
		Points_interface[i].resize(P.Points_interface[i].size());
		for(int j=0; j<P.Points_interface[i].size(); j++ ){
			Points_interface[i][j] = P.Points_interface[i][j];
		}
	}
	
	Points_interface_prev.resize(P.Points_interface_prev.size(), std::vector<Point_3>(0));
	for(int i= 0; i< P.triangles.size(); i++){
		Points_interface_prev[i].resize(P.Points_interface_prev[i].size());
		for(int j=0; j<P.Points_interface_prev[i].size(); j++ ){
			Points_interface_prev[i][j] = P.Points_interface_prev[i][j];
		}
	}
	
	Triangles_interface.resize(P.Triangles_interface.size(), std::vector<Triangle_3>(0));
	for(int i= 0; i< P.triangles.size(); i++){
		Triangles_interface[i].resize(P.Triangles_interface[i].size());
		for(int j=0; j<P.Triangles_interface[i].size(); j++ ){
			Triangles_interface[i][j] = P.Triangles_interface[i][j];
		}
	}
	
	Triangles_interface_prev.resize(P.Triangles_interface_prev.size(), std::vector<Triangle_3>(0));
	for(int i= 0; i< P.triangles.size(); i++){
		Triangles_interface_prev[i].resize(P.Triangles_interface_prev[i].size());
		for(int j=0; j<P.Triangles_interface_prev[i].size(); j++ ){
			Triangles_interface_prev[i][j] = P.Triangles_interface_prev[i][j];
		}
	}
	Position_Triangles_interface.resize(P.Position_Triangles_interface.size(), std::vector< std::vector<int> >(0));
	for(int i= 0; i< P.triangles.size(); i++){
		Position_Triangles_interface[i].resize(P.Position_Triangles_interface[i].size());
		for(int j=0; j<P.Position_Triangles_interface[i].size(); j++ ){
			Position_Triangles_interface[i][j].resize(P.Position_Triangles_interface[i][j].size());
			for(int k=0; k<P.Position_Triangles_interface[i][j].size(); k++ ){
				Position_Triangles_interface[i][j][k] = P.Position_Triangles_interface[i][j][k];
			}
		}
	}
	
	Position_Triangles_interface_prev.resize(P.Position_Triangles_interface_prev.size(), std::vector< std::vector<int> >(0));
	for(int i= 0; i< P.triangles.size(); i++){
		Position_Triangles_interface_prev[i].resize(P.Position_Triangles_interface_prev[i].size());
		for(int j=0; j<P.Position_Triangles_interface_prev[i].size(); j++ ){
			Position_Triangles_interface_prev[i][j].resize(P.Position_Triangles_interface_prev[i][j].size());
			for(int k=0; k<P.Position_Triangles_interface_prev[i][j].size(); k++ ){
				Position_Triangles_interface_prev[i][j][k] = P.Position_Triangles_interface_prev[i][j][k];
			}
		}
	}
	
}
/*!
* \fn void Particule::Affiche()
* \brief Fonction auxiliaire utile pour les tests.
 */
void Particule::Affiche(){
	
	for(int i=0; i<faces.size(); i++){
		cout<<"face "<<i<<endl;
		cout<<" voisin "<<faces[i].voisin<<endl;
		for(int j=0; j<faces[i].size() ;j++){
			cout<<" vertex "<<faces[i].vertex[j].num<<endl;
			for(int k=0; k<faces[i].vertex[j].particules.size();k++){
				cout<<faces[i].vertex[j].particules[k]<<endl;
			}
		}
	}

}

/*!
* \fn void Particule::solve_position(double dt, bool flag_2d)
* \brief Calcul de la position de la particule.
* \param dt pas de temps
* \details  \warning Utilisation de \a Particule.Ff et \a Particule.Mf (forces et moments fluides exerc&eacute;s sur le solide entre t et t+dt/2) \n
  <b>\a Particule.Ff et \a Particule.Mf param&egrave;tres sp&eacute;cifiques au  couplage! </b>
* \return void
 */
void Particule::solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  double eps = 1e-14;//std::numeric_limits<double>::epsilon();
  double rot[3][3];
  if(fixe==1){
    Dx = Vector_3(0.,0.,0.);
    Dxprev = Vector_3(0.,0.,0.);
    u = Vector_3(0.,0.,0.);
    u_half = u;
    e = Vector_3(0.,0.,0.);
    eref = Vector_3(0.,0.,0.);
    rot[0][0] = rot[1][1]= rot[2][2] =1.;
    rot[0][1]=rot[0][2]=rot[1][0]=rot[1][2]=rot[2][0]=rot[2][1]=0.;
    eprev = Vector_3(0.,0.,0.);
    omega = Vector_3(0.,0.,0.);
    omega_half = omega;
  } else {
    if(fixe==0){ //fixe=0: particule mobile
      Dxprev = Dx;
      u = u+(Fi+Ff)/2.*(dt/m);
      u_half = u;
      Dx = Dx+u*dt;
    }
    else if(fixe==2 || fixe==3){//fixe=2: BC en vitesse imposées ! ; fixe=3: fixee en deplacement et rotation seulement selon l'axe y
      //Dx = Vector_3(0.,0.,0.);
      //Dxprev = Vector_3(0.,0.,0.);
      //u = Vector_3(0.,0.,0.);
      //u_half = u;
      Dxprev = Dx;
      //u = u+(Fi+Ff)/2.*(dt/m);
      u_half = u;
      Dx = Dx+u*dt;
    }

    Dx = displacement_BC(x0, Dx, t, T);
      
    //Tests pour verifier qu'on a toujours une matrice de rotation
    for(int i=0;i<3;i++){
      double norm = rotref[i][0]*rotref[i][0]+rotref[i][1]*rotref[i][1]+rotref[i][2]*rotref[i][2];
      if(abs(norm-1.)>eps){
	cout << "Matrice de rotation rotref renomalisee" <<endl;
	cout << "Ligne " << i << " norme=" << norm << endl;
	getchar();
      }
    }
    double vectrot1 = rotref[0][2]-(rotref[1][0]*rotref[2][1]-rotref[2][0]*rotref[1][1]);
    double vectrot2 = rotref[1][2]-(rotref[2][0]*rotref[0][1]-rotref[0][0]*rotref[2][1]);
    double vectrot3 = rotref[2][2]-(rotref[0][0]*rotref[1][1]-rotref[1][0]*rotref[0][1]);
    if(vectrot1*vectrot1+vectrot2*vectrot2+vectrot3*vectrot3>eps){
      cout << "Erreur rotation rotref " << vectrot1 << " " << vectrot2 << " " << vectrot3 << endl;
      //getchar();
    }
    //Calcul de la matrice de rotation totale depuis le repere inertiel jusqu'au temps t et stockage de eprev
    double Q[3][3];
    double eref0 = sqrt(abs(1.-(eref.squared_length())));
    //Recuperation de la matrice de rotation
    Q[0][0] = 1.-2.*(eref[1]*eref[1]+eref[2]*eref[2]);
    Q[0][1] = 2.*(-eref0*eref[2]+eref[0]*eref[1]);
    Q[0][2] = 2.*(eref0*eref[1]+eref[0]*eref[2]);
    Q[1][0] = 2.*(eref0*eref[2]+eref[1]*eref[0]);
    Q[1][1] = 1.-2.*(eref[0]*eref[0]+eref[2]*eref[2]);
    Q[1][2] = 2.*(-eref0*eref[0]+eref[1]*eref[2]);
    Q[2][0] = 2.*(-eref0*eref[1]+eref[2]*eref[0]);
    Q[2][1] = 2.*(eref0*eref[0]+eref[2]*eref[1]);
    Q[2][2] = 1.-2.*(eref[0]*eref[0]+eref[1]*eref[1]);//*/
    double e0 = sqrt(abs(1.-(e.squared_length())));
    /*Recuperation de la matrice de rotation
    rot[0][0] = 1.-2.*(e[1]*e[1]+e[2]*e[2]);
    rot[0][1] = 2.*(-e0*e[2]+e[0]*e[1]);
    rot[0][2] = 2.*(e0*e[1]+e[0]*e[2]);
    rot[1][0] = 2.*(e0*e[2]+e[1]*e[0]);
    rot[1][1] = 1.-2.*(e[0]*e[0]+e[2]*e[2]);
    rot[1][2] = 2.*(-e0*e[0]+e[1]*e[2]);
    rot[2][0] = 2.*(-e0*e[1]+e[2]*e[0]);
    rot[2][1] = 2.*(e0*e[0]+e[2]*e[1]);
    rot[2][2] = 1.-2.*(e[0]*e[0]+e[1]*e[1]);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	Q[i][j] = rot[i][0]*rotref[0][j];
	Q[i][j] += rot[i][1]*rotref[1][j];
	Q[i][j] += rot[i][2]*rotref[2][j];
      }
    }//*/
    /*Impression de la matrice rotref
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	cout << rotref[i][j] << " " ;
      }
      cout << endl;
    }
    getchar();*/
    
    eprev = e;
    //Recuperation du e global a partir de omega
    double Omega[3];
    Omega[0] = Omega[1] = Omega[2] = 0.;
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	       Omega[j] += (omega[k]*Q[k][j]);
      }
    }
    //cout << "debut " << Omega[0] << " " << Omega[1] << " " << Omega[2] << endl;
    //getchar();
    double norm = dt*(abs(Omega[0])+abs(Omega[1])+abs(Omega[2]));
    if(norm>0.25){
      cout << "pas de temps trop grand (solve_position avant moments) : dt=" << dt << " Omega=" << Omega[0] << " " << Omega[1] << " " << Omega[2] << endl;
      getchar();
    }
    double eglob0 = 1.;//sqrt((1+sqrt(1-norm2))/2.);
    double eglob[3];
    for(int j=0;j<3;j++){
      eglob[j] = 0.;//dt*Omega[j]/2./eglob0;
    }
    //Recuperation de la matrice Zn
    double z[3][3];
    z[0][0] = 0.;
    z[0][1] = -Omega[2];
    z[0][2] = Omega[1];
    z[1][0] = Omega[2];
    z[1][1] = 0.;
    z[1][2] = -Omega[0];
    z[2][0] = -Omega[1];
    z[2][1] = Omega[0];
    z[2][2] = 0.;
    //Calcul de la matrice A
    double a[3];
    double d1 = (I[0]+I[1]+I[2])/2.-I[0];
    double d2 = (I[0]+I[1]+I[2])/2.-I[1];
    double d3 = (I[0]+I[1]+I[2])/2.-I[2];
    //Calcul du moment dans le repere inertiel
    double Mx = (Q[0][0]*((Mi+Mf)[0])+Q[1][0]*((Mi+Mf)[1])+Q[2][0]*((Mi+Mf)[2]));
    double My = (Q[0][1]*((Mi+Mf)[0])+Q[1][1]*((Mi+Mf)[1])+Q[2][1]*((Mi+Mf)[2]));
    double Mz = (Q[0][2]*((Mi+Mf)[0])+Q[1][2]*((Mi+Mf)[1])+Q[2][2]*((Mi+Mf)[2]));
    norm = dt*dt/2.*(abs(Mx)/I[0]+abs(My)/I[1]+abs(Mz)/I[2]);
    if(norm>0.25){
      cout << "pas de temps trop grand (solve_position apres moments) : dt=" << dt << " M=" << Mx << " " << My << " " << Mz << " I=" << I[0] << " " << I[1] << " " << I[2] << endl;
      getchar();
    }
    a[0] = -(I[0]*z[1][2]-dt/2.*Mx);
    a[1] = (I[1]*z[0][2]+dt/2.*My);
    a[2] = -(I[2]*z[0][1]-dt/2.*Mz);
    //Resolution du probleme non lineaire
    double etemp0 = 1.;
    double etemp1 = 0.;
    double etemp2 = 0.;
    double etemp3 = 0.;
    double err1 = 1.;
    double err2 = 1.;
    double err3 = 1.;
    double epsilon = 1.e-15;
    int k=0;
    for(k=0; k<1000 && (err1>epsilon*etemp1 || err2>epsilon*etemp2 || err3>epsilon*etemp3); k++){
      if(d2+d3<eps){
	cout << "d2+d3=" << d2+d3 << " I[0]=" << I[0] << endl;
      }
      if(d3+d1<eps){
	cout << "d3+d1=" << d3+d1 << " I[1]=" << I[1] << endl;
      }
      if(d1+d2<eps){
	cout << "d1+d2=" << d1+d2 << " I[2]=" << I[2] << endl;
      }
      double x1 = (dt*a[0]-2.*(d2-d3)*etemp2*etemp3)/(2.*(d2+d3)*etemp0);
      double x2 = (dt*a[1]-2.*(d3-d1)*etemp1*etemp3)/(2.*(d1+d3)*etemp0);
      double x3 = (dt*a[2]-2.*(d1-d2)*etemp1*etemp2)/(2.*(d1+d2)*etemp0);
      if(abs(d2-d3)<epsilon){
	x1 = dt*a[0]/(2.*(d2+d3)*etemp0);
	//cout << "k=" << k << " d2-d3=" << d2-d3 << " I[0]=" << I[0] << endl;
	//getchar();
      }
      if(abs(d3-d1)<epsilon){
	x2 = dt*a[1]/(2.*(d1+d3)*etemp0);
	//cout << "k=" << k << " d3-d1=" << d3-d1 << " I[1]=" << I[1] << endl;
	//getchar();
      }
      if(abs(d1-d2)<epsilon){
	x3 = dt*a[2]/(2.*(d1+d2)*etemp0);
	//cout << "k=" << k << " d1-d2=" << d1-d2 << " I[2]=" << I[2] << endl;
	//getchar();
      }
      etemp1 = x1;
      etemp2 = x2;
      etemp3 = x3;
      if(fixe==3){//fixe=3: on fixe la rotation en x et en z
	etemp1=0.;
	etemp3=0.;
      }
      
// // 	  //Test : on fixe la rotation 
//  			etemp1 = 0.;
//  			etemp2 = 0.;
// 			etemp3 = 0.;
// 			// 			//fin test 
      if(etemp1*etemp1+etemp2*etemp2+etemp3*etemp3>0.5){
	etemp1 /=2.;
	etemp2 /=2.;
	etemp3 /=2.;
	etemp0 = sqrt(1.-etemp1*etemp1-etemp2*etemp2-etemp3*etemp3);
      }
      else{etemp0 = sqrt(1.-etemp1*etemp1-etemp2*etemp2-etemp3*etemp3);}
      err1 = fabs((dt*a[0]-2.*(d2-d3)*etemp2*etemp3)/(2.*(d2+d3)*etemp0)-etemp1);
      err2 = fabs((dt*a[1]-2.*(d3-d1)*etemp1*etemp3)/(2.*(d1+d3)*etemp0)-etemp2);
      err3 = fabs((dt*a[2]-2.*(d1-d2)*etemp1*etemp2)/(2.*(d1+d2)*etemp0)-etemp3);
      if(fixe==3){
	err1 = 0.;
	err3 = 0.;
      }
      
    }//Fin de la boucle de resolution non-lineaire
    if(err1>epsilon || err2>epsilon || err3>epsilon){
      cout << "Probleme de resolution de la rotation, e1=" << etemp1 << " e2=" << etemp2 << " e3=" << etemp3 << endl;
      cout << "erreur=" << err1 << " " << err2 << " " << err3 << endl;
    }
    //cout << k << endl;
    eglob[0] = etemp1;
    eglob[1] = etemp2;
    eglob[2] = etemp3;
    eglob0 = etemp0;
    //Reconstruction de Z^n+1/2
    z[0][0] = (-2.*(eglob[1]*eglob[1]+eglob[2]*eglob[2]))/dt;
    z[0][1] = (-2.*eglob0*eglob[2]+2.*eglob[0]*eglob[1])/dt;
    z[0][2] = (2.*eglob0*eglob[1]+2.*eglob[0]*eglob[2])/dt;
    z[1][0] = (2.*eglob0*eglob[2]+2.*eglob[0]*eglob[1])/dt;
    z[1][1] = (-2.*(eglob[0]*eglob[0]+eglob[2]*eglob[2]))/dt;
    z[1][2] = (-2.*eglob0*eglob[0]+2.*eglob[1]*eglob[2])/dt;
    z[2][0] = (-2.*eglob0*eglob[1]+2.*eglob[0]*eglob[2])/dt;
    z[2][1] = (2.*eglob0*eglob[0]+2.*eglob[1]*eglob[2])/dt;
    z[2][2] = (-2.*(eglob[0]*eglob[0]+eglob[1]*eglob[1]))/dt;
    //Update de la matrice Q
    double Qprev[3][3];
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	       Qprev[i][j] = Q[i][j];
      }
    }
    Q[0][0] = Qprev[0][0]*(1.+dt*z[0][0])+Qprev[0][1]*dt*z[1][0]+Qprev[0][2]*dt*z[2][0];
    Q[0][1] = Qprev[0][0]*dt*z[0][1]+Qprev[0][1]*(1.+dt*z[1][1])+Qprev[0][2]*dt*z[2][1];
    Q[0][2] = Qprev[0][0]*dt*z[0][2]+Qprev[0][1]*dt*z[1][2]+Qprev[0][2]*(1.+dt*z[2][2]);
    Q[1][0] = Qprev[1][0]*(1.+dt*z[0][0])+Qprev[1][1]*dt*z[1][0]+Qprev[1][2]*dt*z[2][0];
    Q[1][1] = Qprev[1][0]*dt*z[0][1]+Qprev[1][1]*(1.+dt*z[1][1])+Qprev[1][2]*dt*z[2][1];
    Q[1][2] = Qprev[1][0]*dt*z[0][2]+Qprev[1][1]*dt*z[1][2]+Qprev[1][2]*(1.+dt*z[2][2]);
    Q[2][0] = Qprev[2][0]*(1.+dt*z[0][0])+Qprev[2][1]*dt*z[1][0]+Qprev[2][2]*dt*z[2][0];
    Q[2][1] = Qprev[2][0]*dt*z[0][1]+Qprev[2][1]*(1.+dt*z[1][1])+Qprev[2][2]*dt*z[2][1];
    Q[2][2] = Qprev[2][0]*dt*z[0][2]+Qprev[2][1]*dt*z[1][2]+Qprev[2][2]*(1.+dt*z[2][2]);
    //Tests pour verifier qu'on a toujours une matrice de rotation
    for(int i=0;i<3;i++){
      double norm = Q[i][0]*Q[i][0]+Q[i][1]*Q[i][1]+Q[i][2]*Q[i][2];
      Q[i][0] /= norm;
      Q[i][1] /= norm;
      Q[i][2] /= norm;
      if(abs(norm-1.)>eps){
	cout << "Matrice de rotation renomalisee" <<endl;
	cout << "Ligne " << i << " norme=" << norm << endl;
	getchar();
      }
    }
    double vect1 = Q[0][2]-(Q[1][0]*Q[2][1]-Q[2][0]*Q[1][1]);
    double vect2 = Q[1][2]-(Q[2][0]*Q[0][1]-Q[0][0]*Q[2][1]);
    double vect3 = Q[2][2]-(Q[0][0]*Q[1][1]-Q[1][0]*Q[0][1]);
    if(vect1*vect1+vect2*vect2+vect3*vect3>eps){
      cout << "Erreur rotation " << vect1 << " " << vect2 << " " << vect3 << endl;
      //getchar();
    }
    //Calcul de eref a partir de Q
    double qref1 = Q[2][1]-Q[1][2];
    double qref2 = Q[0][2]-Q[2][0];
    double qref3 = Q[1][0]-Q[0][1];
    double eref1,eref2,eref3;
    eref1 = signe(qref1)*sqrt(max((1.+Q[0][0]-Q[1][1]-Q[2][2])/4.,0.));
    eref2 = signe(qref2)*sqrt(max((1.+Q[1][1]-Q[0][0]-Q[2][2])/4.,0.));
    eref3 = signe(qref3)*sqrt(max((1.+Q[2][2]-Q[0][0]-Q[1][1])/4.,0.));
    eref = Vector_3(eref1,eref2,eref3);
    /*Version alternative
    Vector_3 qref = Vector_3(Q[2][1]-Q[1][2], Q[0][2]-Q[2][0], Q[1][0]-Q[0][1])/4.;
    eref0 = sqrt((1+sqrt(1-4.*qref.squared_length()))/2.);
    eref = qref/eref0;//*/
    //Recuperation de la matrice de rotation de la particule
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	rot[i][j] = Q[i][0]*rotref[j][0];
	rot[i][j] += Q[i][1]*rotref[j][1];
	rot[i][j] += Q[i][2]*rotref[j][2];
      }
    }
    //Calcul de e a partir de rot
    double q1 = rot[2][1]-rot[1][2];
    double q2 = rot[0][2]-rot[2][0];
    double q3 = rot[1][0]-rot[0][1];
    double e1,e2,e3;
    e1 = signe(q1)*sqrt(max((1.+rot[0][0]-rot[1][1]-rot[2][2])/4.,0.));
    e2 = signe(q2)*sqrt(max((1.+rot[1][1]-rot[0][0]-rot[2][2])/4.,0.));
    e3 = signe(q3)*sqrt(max((1.+rot[2][2]-rot[0][0]-rot[1][1])/4.,0.));
    e = Vector_3(e1,e2,e3);
    //*/
    /*Version alternative
    Vector_3 q = Vector_3(rot[2][1]-rot[1][2], rot[0][2]-rot[2][0], rot[1][0]-rot[0][1])/4.;
    e0 = sqrt((1+sqrt(1-4.*q.squared_length()))/2.);
    e = q/e0;//*/
    
    //Calcul de Omega^n+1/2
    double omega1 = 0.;
    for(int i=0;i<3;i++){
      //for(int j=0;j<3;j++){
      //omega1 -= 1./2.*Qprev[1][i]*z[i][j]*(Qprev[2][j]+Q[2][j]);
      omega1 += 1./2./dt*(Q[2][i]*Qprev[1][i]-Qprev[2][i]*Q[1][i]);
      //}
    }
    double omega2 = 0.;
    for(int i=0;i<3;i++){
      //for(int j=0;j<3;j++){
      //omega2 += 1./2.*Qprev[0][i]*z[i][j]*(Qprev[2][j]+Q[2][j]);
      omega2 += 1./2./dt*(Q[0][i]*Qprev[2][i]-Qprev[0][i]*Q[2][i]);
      //}
    }
    double omega3 = 0.;
    for(int i=0;i<3;i++){
      //for(int j=0;j<3;j++){
      //omega3 -= 1./2.*Qprev[0][i]*z[i][j]*(Qprev[1][j]+Q[1][j]);
      omega3 += 1./2./dt*(Q[1][i]*Qprev[0][i]-Qprev[1][i]*Q[0][i]);
      //}
    }
    omega = Vector_3(omega1,omega2,omega3);
    omega_half = omega;
		if(flag_2d){
			omega = Vector_3(0.,0.,omega3);
			omega_half = omega;
		}
  }//Fin du calcul dans le cas d'une particule libre
  /*Test de fixer la rotation
  rot[0][0]= rot[1][1] = rot[2][2] =1.;
  rot[0][1] = rot[0][2] =rot[1][0] = rot[1][2] = rot[2][0] = rot[2][1] = 0.;
  e = Vector_3(0,0,0);
  //rotprev[0][0]= rotprev[1][1] = rotprev[2][2] =1.;
  //rotprev[0][1] = rotprev[0][2] =rotprev[1][0] = rotprev[1][2] = rotprev[2][0] = rotprev[2][1] = 0.;

  omega = Vector_3(0.,0.,0.);
  omega_half = omega;
  //fin test */

  //Mise a jour de la transformation donnant le mouvement de la particule
  mvt_tprev = mvt_t;
  Aff_transformation_3 rotation(rot[0][0],rot[1][0],rot[2][0],rot[0][1],rot[1][1],rot[2][1],rot[0][2],rot[1][2],rot[2][2]);
  Aff_transformation_3 translation(Vector_3(Point_3(0.,0.,0.),x0)+Dx);
  Aff_transformation_3 translation_inv(Vector_3(x0,Point_3(0.,0.,0.)));
  mvt_t = translation*(rotation*translation_inv);
	//cout<<"position du centre de la particule "<<x0+Dx<<endl;
}

/*!
* \fn void Particule::solve_vitesse(double dt, bool flag_2d)
* \brief Calcul de la vitesse de la particule.
* \param dt pas de temps
\details  \warning Utilisation de \a Particule.Ff et \a Particule.Mf (forces et moments fluides exerc&eacute;s sur le solide entre t et t+dt/2) \n
<b>\a Particule.Ff et \a Particule.Mf param&egrave;tres sp&eacute;cifiques au  couplage! </b>
* \return void
 */
void Particule::solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  if(fixe==1){
    u = Vector_3(0.,0.,0.);
    omega = Vector_3(0.,0.,0.);
  } else {
    if(fixe==0){
      u = u+(Fi+Ff)/2.*(dt/m)*Amort;// + velocity_BC(x0, t, T, Dx); //Conditions aux limites en vitesse ajoutées ici
    }
    else if(fixe==2 || fixe==3){
      //u = velocity_BC(x0, t, T, Dx);
      u = u+(Fi+Ff)/2.*(dt/m);
      u.vec[2] = velocity_BC_bis(x0, t, T, Dx);
    }
    
    //Calcul de la matrice de rotation totale depuis le repï¿½re inertiel jusqu'au temps t
    double Q[3][3];
    double eref0 = sqrt(abs(1.-(eref.squared_length())));
    //Recuperation de la matrice de rotation
    Q[0][0] = 1.-2.*(eref[1]*eref[1]+eref[2]*eref[2]);
    Q[0][1] = 2.*(-eref0*eref[2]+eref[0]*eref[1]);
    Q[0][2] = 2.*(eref0*eref[1]+eref[0]*eref[2]);
    Q[1][0] = 2.*(eref0*eref[2]+eref[1]*eref[0]);
    Q[1][1] = 1.-2.*(eref[0]*eref[0]+eref[2]*eref[2]);
    Q[1][2] = 2.*(-eref0*eref[0]+eref[1]*eref[2]);
    Q[2][0] = 2.*(-eref0*eref[1]+eref[2]*eref[0]);
    Q[2][1] = 2.*(eref0*eref[0]+eref[2]*eref[1]);
    Q[2][2] = 1.-2.*(eref[0]*eref[0]+eref[1]*eref[1]);//*/
    double e0 = sqrt(abs(1.-(e.squared_length())));
    /*Recuperation de la matrice de rotation
    double rot[3][3];
    rot[0][0] = 1.-2.*(e[1]*e[1]+e[2]*e[2]);
    rot[0][1] = 2.*(-e0*e[2]+e[0]*e[1]);
    rot[0][2] = 2.*(e0*e[1]+e[0]*e[2]);
    rot[1][0] = 2.*(e0*e[2]+e[1]*e[0]);
    rot[1][1] = 1.-2.*(e[0]*e[0]+e[2]*e[2]);
    rot[1][2] = 2.*(-e0*e[0]+e[1]*e[2]);
    rot[2][0] = 2.*(-e0*e[1]+e[2]*e[0]);
    rot[2][1] = 2.*(e0*e[0]+e[2]*e[1]);
    rot[2][2] = 1.-2.*(e[0]*e[0]+e[1]*e[1]);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	Q[i][j] = rot[i][0]*rotref[0][j];
	Q[i][j] += rot[i][1]*rotref[1][j];
	Q[i][j] += rot[i][2]*rotref[2][j];
      }
    } //*/
    //Recuperation de Zn+1/2 a partir de omega
    double Omega[3];
    Omega[0] = Omega[1] = Omega[2] = 0.;
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	Omega[j] += (omega[k]*Q[k][j]);
      }
    }
    //cout << "debut " << Omega[0] << " " << Omega[1] << " " << Omega[2] << endl;
    //getchar();
    double norm2 = dt*dt*(Omega[0]*Omega[0]+Omega[1]*Omega[1]+Omega[2]*Omega[2]);
    if(norm2>1.){
      cout << "pas de temps trop grand (solve_vitesse) : dt=" << dt << " Omega=" << Omega[0] << " " << Omega[1] << " " << Omega[2] << endl;
      getchar();
    }
    double eglob0 = sqrt((1.+sqrt(1.-norm2))/2.);
    double eglob[3];
    for(int j=0;j<3;j++){
      eglob[j] = dt*Omega[j]/2./eglob0;
    }
    //Recuperation de la matrice Zn+1/2
    //double eglob0 = sqrt(1.-eglob[0]*eglob[0]-eglob[1]*eglob[1]-eglob[2]*eglob[2]);
    double z[3][3];
    z[0][0] = (-2.*(eglob[1]*eglob[1]+eglob[2]*eglob[2]))/dt;
    z[0][1] = (-2.*eglob0*eglob[2]+2.*eglob[0]*eglob[1])/dt;
    z[0][2] = (2.*eglob0*eglob[1]+2.*eglob[0]*eglob[2])/dt;
    z[1][0] = (2.*eglob0*eglob[2]+2.*eglob[0]*eglob[1])/dt;
    z[1][1] = (-2.*(eglob[0]*eglob[0]+eglob[2]*eglob[2]))/dt;
    z[1][2] = (-2.*eglob0*eglob[0]+2.*eglob[1]*eglob[2])/dt;
    z[2][0] = (-2.*eglob0*eglob[1]+2.*eglob[0]*eglob[2])/dt;
    z[2][1] = (2.*eglob0*eglob[0]+2.*eglob[1]*eglob[2])/dt;
    z[2][2] = (-2.*(eglob[0]*eglob[0]+eglob[1]*eglob[1]))/dt;
    //Calcul de la matrice A
    double a[3];
    double d1 = (I[0]+I[1]+I[2])/2.-I[0];
    double d2 = (I[0]+I[1]+I[2])/2.-I[1];
    double d3 = (I[0]+I[1]+I[2])/2.-I[2];
    //Calcul du moment dans le repere inertiel
    double Mx = (Q[0][0]*((Mi+Mf)[0])+Q[1][0]*((Mi+Mf)[1])+Q[2][0]*((Mi+Mf)[2]));
    double My = (Q[0][1]*((Mi+Mf)[0])+Q[1][1]*((Mi+Mf)[1])+Q[2][1]*((Mi+Mf)[2]));
    double Mz = (Q[0][2]*((Mi+Mf)[0])+Q[1][2]*((Mi+Mf)[1])+Q[2][2]*((Mi+Mf)[2]));
    a[0] = -(d2*z[1][2]-d3*z[2][1]-dt/2.*Mx);
    a[1] = (d1*z[0][2]-d3*z[2][0]+dt/2.*My);
    a[2] = -(d1*z[0][1]-d2*z[1][0]-dt/2.*Mz);
    //Resolution du probleme lineaire sur Zn+1
    z[0][0] = 0.;
    z[0][1] = -a[2]/I[2];
    z[0][2] = a[1]/I[1];
    z[1][0] = -z[0][1];
    z[1][1] = 0.;
    z[1][2] = -a[0]/I[0];
    z[2][0] = -z[0][2];
    z[2][1] = -z[1][2];
    z[2][2]= 0.;
    //Calcul de Omega^n+1
    double omega1 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega1 -= Q[1][i]*z[i][j]*Q[2][j];
      }
    }
    double omega2 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega2 += Q[0][i]*z[i][j]*Q[2][j];
      }
    }
    double omega3 = 0.;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	omega3 -= Q[0][i]*z[i][j]*Q[1][j];
      }
    }
//     //Test pour fixer les composantes x et y de la rotation
			if(flag_2d){
			 		omega1 = 0.;
			 		omega2 = 0.;
			}
			if(fixe==3){
			  omega1=0.;
			  omega3=0.;
			}
			
// 		//fin test 
    omega = Vector_3(omega1,omega2,omega3);
    /*Test de fixer la rotation
		rot[0][0]= rot[1][1] = rot[2][2] =1.;
		rot[0][1] = rot[0][2] =rot[1][0] = rot[1][2] = rot[2][0] = rot[2][1] = 0.;
		e = Vector_3(0,0,0);
		//rotprev[0][0]= rotprev[1][1] = rotprev[2][2] =1.;
		//rotprev[0][1] = rotprev[0][2] =rotprev[1][0] = rotprev[1][2] = rotprev[2][0] = rotprev[2][1] = 0.;
		omega = Vector_3(0.,0.,0.);
		omega_half = omega;
		//fin test */
  }//Fin du calcul dans le cas d'une particule libre
}

/*void Particule::solve_vitesse_plas(const double& dt, const bool& flag_2d){
  if(fixe==1){
    u_plas = Vector_3(0.,0.,0.);
  } else {
    if(fixe==0){
      u_plas = u_plas + Fi_plas/2.*(dt/m);
    }
    else if(fixe==2 || fixe==3){
      u_plas = Vector_3(0.,0.,0.);
    }
  }
  }*/

/*!
* \fn double Particule::volume()
* \brief Fonction auxilaire utile pour les tests. Calcul du volume de la particule.
* \return double
 */
// double Particule::volume(){
// 	
//   double vol = 0.;
//   std::vector<Point_3> Points_poly; 
// 	
//   for(int l= 0; l<triangles.size(); l++)
//   {
//     Points_poly.push_back(triangles[l][0]);
//     Points_poly.push_back(triangles[l][1]);
//     Points_poly.push_back(triangles[l][2]);
//   }	
//   Finite_cells_iterator cit;
//   Triangulation T(Points_poly.begin(), Points_poly.end());
// 	
//   for (cit = T.finite_cells_begin(); cit!= T.finite_cells_end(); cit++){
//     vol+= CGAL::to_double(T.tetrahedron( cit).volume());
//   }
// 	
//   return vol;
// }
double Particule::volume(){
	
	double vol = 0.;
	
	Point_3 center;
	
	std::vector<Point_3> Points_poly; 
	
	for(int l= 0; l<triangles.size(); l++)
	{
		Points_poly.push_back(triangles[l][0]);
		Points_poly.push_back(triangles[l][1]);
		Points_poly.push_back(triangles[l][2]);
	}	
	center = centroid(Points_poly.begin(),Points_poly.end());
	
	for(int l= 0; l<triangles.size(); l++)
	{
		Tetrahedron tetra(center, triangles[l][0], triangles[l][1], triangles[l][2]);
		vol += abs((tetra.volume()));
	}
	
	return vol;
}

/*!
* \fn Particule::vitesse_parois(Point_3& X_f)
* \brief Vitesse au centre de la parois au temps t. \n
* \f$ V_f = V_I + \Omega_{rot} \wedge \left( X_f - X_I \right). \f$ 
* \f$ V_I \f$ -vitesse de la particule(\a Particule.u_half),  \f$ X_I \f$ -centre de la particule(\a Particule.x0 + \a Particule.Dx),  \f$ \Omega_{rot} \f$ -rotation de la particule(\a Particule.omega_half).
* \param X_f centre de la parois
* \warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
* \return Vector_3
*/
Vector_3 Particule::vitesse_parois(const Point_3& X_f){
		
  Vector_3 V_f = u_half + cross_product(omega_half, Vector_3(x0 + Dx,X_f));

	return V_f;
}	
/*!
* \fn Particule::vitesse_parois_prev(Point_3& X_f)
* \brief Vitesse au centre de la parois au temps t-dt.
 * \f$ V_f = V_I + \Omega_{rot} \wedge \left( X_f - X_I \right). \f$ \n
 * \f$ V_I \f$ -vitesse de la particule(\a Particule.u_half),  \f$ X_I \f$ -centre de la particule(\a Particule.x0 + \a Particule.Dxprev),  \f$ \Omega_{rot} \f$ -rotation de la particule(\a Particule.omega_half).
* \param X_f centre de la parois
* \warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
* \return Vector_3
 */
Vector_3 Particule::vitesse_parois_prev(const Point_3& X_f){
	
  Vector_3 V_f = u_half + cross_product(omega_half, Vector_3(x0 + Dxprev,X_f));
	
	return V_f;
}	

/*!
* \fn void Face::compProjectionIntegrals(double &P1, double &Pa, double &Pb, double &Paa, double &Pab, double &Pbb, double &Paaa, double &Paab, double &Pabb, double &Pbbb, int a, int b, int c)
*\brief Calcul des projections.
* \details Utilisation de la fonction d&eacute;crite par Brian Mirtich 1996(cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar).
* \warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
* \return void
*/
void Face::compProjectionIntegrals(double &P1, double &Pa, double &Pb, double &Paa, double &Pab, double &Pbb, double &Paaa, double &Paab, double &Pabb, double &Pbbb, const int& a, const int& b, const int& c){
  //Utilisation de la fonction decrite par Brian Mirtich 1996 (cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar)
  P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.;
  for(int i=0;i<size();i++){
    int j= (i+1)%(size());
    double a0 = (vertex[i].pos[a]);
    double b0 = (vertex[i].pos[b]);
    double a1 = (vertex[j].pos[a]);
    double b1 = (vertex[j].pos[b]);
    double Da = a1-a0;
    double Db = b1-b0;
    double a02 = a0*a0;
    double a03 = a0*a02;
    double a04 = a0*a03;
    double b02 = b0*b0;
    double b03 = b0*b02;
    double b04 = b0*b03;
    double a12 = a1*a1;
    double a13 = a1*a12;
    double b12 = b1*b1;
    double b13 = b1*b12;
    double C1 = a1+a0;
    double Ca = a1*C1+a02;
    double Caa = a1*Ca+a03;
    double Caaa = a1*Caa+a04;
    double Cb = b12+b1*b0+b02;
    double Cbb = b1*Cb+b03;
    double Cbbb = b1*Cbb+b04;
    double Cab = 3.*a12+2.*a1*a0+a02;
    double Kab = a12+2.*a1*a0+3.*a02;
    double Caab = a0*Cab+4.*a13;
    double Kaab = a1*Kab+4.*a03;
    double Cabb = 4.*b13+3.*b12*b0+2.*b1*b02+b03;
    double Kabb = b13+2.*b12*b0+3.*b1*b02+4.*b03;
    P1 += Db*C1;
    Pa += Db*Ca; Paa += Db*Caa; Paaa += Db*Caaa;
    Pb += Da*Cb; Pbb += Da*Cbb; Pbbb += Da*Cbbb;
    Pab += Db*(b1*Cab+b0*Kab);
    Paab += Db*(b1*Caab+b0*Kaab);
    Pabb += Da*(a1*Cabb+a0*Kabb);
  }
  P1 /= 2.;
  Pa /= 6.; Paa /= 12.; Paaa /= 20.;
  Pb /=-6.; Pbb /=-12.; Pbbb /=-20.;
  Pab /= 24.;
  Paab /= 60.;
  Pabb /= -60.;
}
/*!
* \fn void Face::compFaceIntegrals(double &Fa, double &Fb, double &Fc, double &Faa, double &Fbb, double &Fcc, double &Faaa, double &Fbbb, double &Fccc, double &Faab, double &Fbbc, double &Fcca, double na, double nb, double nc, int a, int b, int c)
* \brief Calcul des int&eacute;grales sur les faces. 
* \details Utilisation de la fonction d&eacute;crite par Brian Mirtich 1996(cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar).
* \warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
* \return void
*/
void Face::compFaceIntegrals(double &Fa, double &Fb, double &Fc, double &Faa, double &Fbb, double &Fcc, double &Faaa, double &Fbbb, double &Fccc, double &Faab, double &Fbbc, double &Fcca, const double& na, const double& nb, const double& nc, const int& a, const int& b, const int& c){
  //Utilisation de la fonction decrite par Brian Mirtich 1996 (cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar)
  double P1,Pa,Pb,Paa,Pab,Pbb,Paaa,Paab,Pabb,Pbbb;
  Vector_3 p(Point_3(0,0,0),vertex[0].pos);
  double w = -(normale*p);
  double k1 = 1./nc;
  double k2 = k1*k1;
  double k3 = k1*k2;
  double k4 = k1*k3;
  compProjectionIntegrals(P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb, a, b, c);
  Fa = k1*Pa;
  Fb = k1*Pb;
  Fc = -k2*(na*Pa+nb*Pb+w*P1);
  Faa = k1*Paa;
  Fbb = k1*Pbb;
  Fcc = k3*(na*na*Paa+2.*na*nb*Pab+nb*nb*Pbb+2.*na*w*Pa+2.*nb*w*Pb+w*w*P1);
  Faaa = k1*Paaa;
  Fbbb = k1*Pbbb;
  Fccc = -k4*(na*na*na*Paaa+3.*na*na*nb*Paab+3.*na*nb*nb*Pabb+nb*nb*nb*Pbbb+3.*na*na*w*Paa+6.*na*nb*w*Pab+3.*nb*nb*w*Pbb+3.*na*w*w*Pa+3.*nb*w*w*Pb+w*w*w*P1);
  Faab = k1*Paab;
  Fbbc = -k2*(na*Pabb+nb*Pbbb+w*Pbb);
  Fcca = k3*(na*na*Paaa+2.*na*nb*Paab+nb*nb*Pabb+2.*na*w*Paa+2.*nb*w*Pab+w*w*Pa);
}

/*!
* \fn void Particule::CompVolumeIntegrals(double &T1, double &Tx, double &Ty, double &Tz, double &Txx, double &Tyy, double &Tzz, double &Txy, double &Tyz, double &Tzx)
* \brief Calcul des int&eacute;grales de volume.
*\details Utilisation de la fonction d&eacute;crite par Brian Mirtich 1996(cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar).
* \warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
* \return void
*/
void Particule::CompVolumeIntegrals(double &T1, double &Tx, double &Ty, double &Tz, double &Txx, double &Tyy, double &Tzz, double &Txy, double &Tyz, double &Tzx){
  //Utilisation de la fonction decrite par Brian Mirtich 1996 (cf www.cs.berkeley.edu/~jfc/mirtich/code/volumeIntegration.tar)
  T1 = Tx=Ty=Tz=Txx=Tyy=Tzz=Txy=Tyz=Tzx=0.;
  for(int i=0;i<faces.size();i++){
    double Fx,Fy,Fz,Fxx,Fyy,Fzz,Fxxx,Fyyy,Fzzz,Fxxy,Fyyz,Fzzx;
    double nx,ny,nz,na,nb,nc;
    nx=(faces[i].normale[0]);
    ny=(faces[i].normale[1]);
    nz=(faces[i].normale[2]);
    //Choix d'une permutation orientee abc telle que nc soit maximale
    if(abs(nx)>abs(ny)){
      if(abs(nx)>abs(nz)){
	//Cas a=y,b=z,c=x
	na = ny; nb = nz; nc = nx;
	faces[i].compFaceIntegrals(Fy, Fz, Fx, Fyy, Fzz, Fxx, Fyyy, Fzzz, Fxxx, Fyyz, Fzzx, Fxxy,na,nb,nc,1,2,0);
      }
      else {
	//Cas a=x,b=y,c=z
	na = nx; nb = ny; nc = nz;
	faces[i].compFaceIntegrals(Fx, Fy, Fz, Fxx, Fyy, Fzz, Fxxx, Fyyy, Fzzz, Fxxy, Fyyz, Fzzx,na,nb,nc,0,1,2);
      }
    }
    else {
      if(abs(ny)>abs(nz)){
	//Cas a=z,b=x,c=y
	na = nz; nb = nx; nc = ny;
	faces[i].compFaceIntegrals(Fz, Fx, Fy, Fzz, Fxx, Fyy, Fzzz, Fxxx, Fyyy, Fzzx, Fxxy, Fyyz,na,nb,nc,2,0,1);
      }
      else{
	//Cas a=x,b=y,c=z
	na = nx; nb = ny; nc = nz;
	faces[i].compFaceIntegrals(Fx, Fy, Fz, Fxx, Fyy, Fzz, Fxxx, Fyyy, Fzzz, Fxxy, Fyyz, Fzzx,na,nb,nc,0,1,2);
      }
    }
    //Calcul des integrales
    T1 += nx*Fx;
    Tx += nx*Fxx;
    Ty += ny*Fyy;
    Tz += nz*Fzz;
    Txx += nx*Fxxx;
    Tyy += ny*Fyyy;
    Tzz += nz*Fzzz;
    Txy += nx*Fxxy;
    Tyz += ny*Fyyz;
    Tzx += nz*Fzzx;
  }
  Tx /= 2.;
  Ty /= 2.;
  Tz /= 2.;
  Txx /=3.;
  Tyy /=3.;
  Tzz /=3.;
  Txy /=2.;
  Tyz /=2.;
  Tzx /=2.;
}

struct Mat3x3
{
  double tab[3][3];
};

struct Vect3
{
  double vec[3];
};

//Fonction rot pour la routine jacobi3x3
inline void rot(Mat3x3 &a, const double& s, const double& tau, const int& i, const int& j, const int& k, const int& l)
{
  double g,h;
  
  g = a.tab[i][j];
  h = a.tab[k][l];
  a.tab[i][j] = g-s*(h+g*tau);
  a.tab[k][l] = h+s*(g-h*tau);
}

//Diagonalisation de la matrice 3x3 a par la methode de Jacobi
//cf Numerical Recipes C++
//a est la matrice qu'on diagonalise, d la diagonale des valeurs propres, v la matrice des vecteurs propres, nrot le nombre d'iterations de Jacobi
void jacobi3x3(Mat3x3 &a, Vect3 &d, Mat3x3 &v, int &nrot)
{
  int i,j,ip,iq;
  double tresh,theta,tau,t,sm,s,h,g,c;
  
  const int n=3;
  double b[n],z[n];
  //Initialisation de v a l'identite
  for(ip=0;ip<n;ip++){
    for(iq=0;iq<n;iq++){
      v.tab[ip][iq]=0.;
    }
    v.tab[ip][ip]=1.;
  }
  //Initialisation de b et d a la diagonale de a
  for(ip=0;ip<n;ip++){
    b[ip]=d.vec[ip]=a.tab[ip][ip];
    z[ip]=0.;
  }
  nrot = 0;
  for(i=1;i<=50;i++){
    sm=0.;
    //Somme des magnitudes des elements hors diagonale
    for(ip=0;ip<n-1;ip++){
      for(iq=ip+1;iq<n;iq++){
	sm+=fabs(a.tab[ip][iq]);
      }
    }
    //Si on a convergence exacte
    if(sm==0.){
      return;
    }
    //Sur les trois premiers sweeps
    if(i<4){
      tresh=0.2*sm/(n*n);
    } 
    //et ensuite...
    else {
      tresh=0.;
    }
    //Debut du sweep
    for(ip=0;ip<n-1;ip++){
      for(iq=ip+1;iq<n;iq++){
	g = 100.*fabs(a.tab[ip][iq]);
	//Apres 4 sweeps, sauter la rotation si l'element off-diagonal est petit
	if(i>4 && (fabs(d.vec[ip])+g)==fabs(d.vec[ip]) && (fabs(d.vec[iq])+g)==fabs(d.vec[iq])){
	  a.tab[ip][iq]=0.;
	}
	else if(fabs(a.tab[ip][iq])>tresh){
	  h=d.vec[iq]-d.vec[ip];
	  if((fabs(h)+g)==fabs(h)){
	    t=(a.tab[ip][iq])/h;
	  } else {
	    theta=0.5*h/(a.tab[ip][iq]);
	    t = 1./(fabs(theta)+sqrt(1.+theta*theta));
	    if(theta<0.){
	      t = -t;
	    }
	  }
	  c = 1./sqrt(1+t*t);
	  s = t*c;
	  tau = s/(1.+c);
	  h = t*a.tab[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d.vec[ip] -= h;
	  d.vec[iq] += h;
	  a.tab[ip][iq] = 0.;
	  for(j=0;j<ip;j++){
	    rot(a,s,tau,j,ip,j,iq);
	  }
	  for(j=ip+1;j<iq;j++){
	    rot(a,s,tau,ip,j,j,iq);
	  }
	  for(j=iq+1;j<n;j++){
	    rot(a,s,tau,ip,j,iq,j);
	  }
	  for(j=0;j<n;j++){
	    rot(v,s,tau,j,ip,j,iq);
	  }
	  ++nrot;
	}
      }
    }//Fin du sweep
    for(ip=0;ip<n;ip++){
      b[ip] += z[ip];
      d.vec[ip] = b[ip];
      z[ip] = 0.;
    }
  }
  cout << "Trop grand nombre d'iterations de la routine jacobi3x3" << endl;
}


/*!
* \fn void Particule::Inertie(const double &rho)
* \brief Calcul d'inertie de la particule. 
* \warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
* \return void
*/
void Particule::Inertie(const double &rho){
  double eps = 1.e-14;//std::numeric_limits<double>::epsilon();
  double T1,Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Tzx;
  CompVolumeIntegrals(T1,Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Tzx);
  double R[3][3];
  double xG = (x0[0]);
  double yG = (x0[1]);
  double zG = (x0[2]);
  
//  cout << "T=" << endl;
//  cout << Txx << " " << Txy << " " << Tzx << endl;
//  cout << Txy << " " << Tyy << " " << Tyz << endl;
//  cout << Tzx << " " << Tyz << " " << Tzz << endl;
//  getchar();
  
  R[0][0] = rho*(Tyy-2.*yG*Ty+yG*yG*T1+Tzz-2.*zG*Tz+zG*zG*T1);
  R[1][0] = R[0][1] = rho*(Txy-yG*Tx-xG*Ty+xG*yG*T1);
  R[2][0] = R[0][2] = rho*(Tzx-zG*Tx-xG*Tz+xG*zG*T1);
  R[1][1] = rho*(Txx-2.*xG*Tx+xG*xG*T1+Tzz-2.*zG*Tz+zG*zG*T1);
  R[1][2] = R[2][1] = rho*(Tyz-zG*Ty-yG*Tz+yG*zG*T1);
  R[2][2] = rho*(Tyy-2.*yG*Ty+yG*yG*T1+Txx-2.*xG*Tx+xG*xG*T1);
  
//  cout << "R=" << endl;
//  cout << R[0][0] << " " << R[0][1] << " " << R[0][2] << endl;
//  cout << R[1][0] << " " << R[1][1] << " " << R[1][2] << endl;
//  cout << R[2][0] << " " << R[2][1] << " " << R[2][2] << endl;
//  getchar();

  //Masse et volume
  V = T1;
  m = rho*T1;
  if(m<eps){
    cout<< "masse nulle " << m << endl;
    getchar();
  }
  //Calcul des moments d'inertie
  //Nouvelle version utilisant la methode de Jacobi (Numerical Recipes c++)
  //if(max(max(abs(R[0][1]),abs(R[0][2])),abs(R[1][2]))>eps){
    Mat3x3 A,V;
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	A.tab[i][j] = R[i][j];
      }
    }
    Vect3 d;
    int n=0;
    jacobi3x3(A,d,V,n);
    for(int i=0;i<3;i++){
      I[i] = d.vec[i];
      for(int j=0;j<3;j++){
	rotref[i][j] = V.tab[j][i];
      }
    }
    /*} else {
    I[0] = R[0][0];
    I[1] = R[1][1];
    I[2] = R[2][2];
    rotref[0][0] =rotref[1][1] = rotref[2][2] = 1.;
    rotref[0][1] = rotref[0][2] = rotref[1][0] = rotref[1][2] = rotref[2][0] = rotref[2][1] = 0.;
    }*/
    //  cout << "I=" << I[0] << " " << I[1] << " " << I[2] << endl;
//  getchar();
  
  //Test : produit scalaire des deux premieres colonnes
  double scal = rotref[0][0]*rotref[0][1]+rotref[1][0]*rotref[1][1]+rotref[2][0]*rotref[2][1];
// 			cout << "produit scalaire " << scal << endl;
  if(abs(scal)>eps){
    cout << "scal=" << scal << endl;
  }
  
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      if(rotref[i][j]!=rotref[i][j]){
	cout << "rotref "<< rotref[i][j] << " " << i << " " << j << endl;
	getchar();
      }
    }
  }
  
  //Test sur le determinant de la matrice de rotation (1 ou -1)
  double det = rotref[0][2]*(rotref[1][0]*rotref[2][1]-rotref[2][0]*rotref[1][1]);
  det += rotref[1][2]*(rotref[2][0]*rotref[0][1]-rotref[0][0]*rotref[2][1]);
  det += rotref[2][2]*(rotref[0][0]*rotref[1][1]-rotref[1][0]*rotref[0][1]);
  if(det<0.){
    for(int i=0;i<3;i++){
      rotref[i][2] *= -1.;
    }
  }
  for(int i=0;i<3;i++){
    int j = (i+1)%3;
    if(abs(rotref[0][i]*rotref[0][j]+rotref[1][i]*rotref[1][j]+rotref[2][i]*rotref[2][j])>eps){
      cout << "erreur dans le calcul des moments d'inertie" << endl;
      getchar();
    }
  }

  //Calcul de eref dans le référentiel de la particule
  double Q[3][3];
  double e0 = sqrt(abs(1.-(e.squared_length())));
  double rot[3][3];
  //Recuperation de la matrice de rotation
  rot[0][0] = 1.-2.*(e[1]*e[1]+e[2]*e[2]);
  rot[0][1] = 2.*(-e0*e[2]+e[0]*e[1]);
  rot[0][2] = 2.*(e0*e[1]+e[0]*e[2]);
  rot[1][0] = 2.*(e0*e[2]+e[1]*e[0]);
  rot[1][1] = 1.-2.*(e[0]*e[0]+e[2]*e[2]);
  rot[1][2] = 2.*(-e0*e[0]+e[1]*e[2]);
  rot[2][0] = 2.*(-e0*e[1]+e[2]*e[0]);
  rot[2][1] = 2.*(e0*e[0]+e[2]*e[1]);
  rot[2][2] = 1.-2.*(e[0]*e[0]+e[1]*e[1]);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Q[i][j] = rot[i][0]*rotref[0][j];
      Q[i][j] += rot[i][1]*rotref[1][j];
      Q[i][j] += rot[i][2]*rotref[2][j];
    }
  }
  //Calcul de eref a partir de Q
  double qref1 = Q[2][1]-Q[1][2];
  double qref2 = Q[0][2]-Q[2][0];
  double qref3 = Q[1][0]-Q[0][1];
  double eref1,eref2,eref3;
  eref1 = signe(qref1)*sqrt(max((1.+Q[0][0]-Q[1][1]-Q[2][2])/4.,0.));
  eref2 = signe(qref2)*sqrt(max((1.+Q[1][1]-Q[0][0]-Q[2][2])/4.,0.));
  eref3 = signe(qref3)*sqrt(max((1.+Q[2][2]-Q[0][0]-Q[1][1])/4.,0.));
  eref = Vector_3(eref1,eref2,eref3);
  /*Version alternative
  Vector_3 qref = Vector_3(Q[2][1]-Q[1][2], Q[0][2]-Q[2][0], Q[1][0]-Q[0][1])/4.;
  double eref0 = sqrt((1+sqrt(1-4.*qref.squared_length()))/2.);
  eref = qref/eref0;//*/
  
  //Calcul des moments d'inertie des faces (pour le calcul des torsions)
  for(int i=0;i<faces.size();i++){
    faces[i].Inertie();
  }
}
/*!
* \fn void Particule::Volume_libre()
* \brief Calcul du volume libre. 
* \warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
* \return void
*/
void Particule::Volume_libre(){
  Vl = 0.;
  for(int i=0;i<faces.size();i++){
    if(faces[i].voisin == -1){
      Vector_3 v1(faces[i].vertex[0].pos,faces[i].vertex[1].pos);
      Vector_3 v2(faces[i].vertex[0].pos,faces[i].vertex[2].pos);
      Vector_3 v3(x0,faces[i].vertex[0].pos);
      Vl += 1./6.*(cross_product(v1,v2)*v3);
    }
  }
}

/*!
*\fn Solide::Solide()
*\brief Constructeur par d&eacute;faut. 
*/
Solide::Solide(const double& E, const double& nu){
  lambda = E * nu / (1+nu) / (1 - 2.*nu);
  mu = E / 2. / (1+nu);
}

Solide::Solide(){
  lambda = 0.;
  mu = 0.;
}

Solide::Solide(const std::vector<Particule> & Part){
  lambda = 0.;
  mu = 0.;
  for(int i=0; i<Part.size(); i++){
    solide.push_back(Part[i]);
  }
}
/*!
*\fn Cellule::~Cellule()
*\brief Destructeur.
*/ 
Solide::~Solide(){   
}

/*!
*\fn Solide & Solide:: operator=(const Solide &S)
*\brief op&eacute;rateur = Surcharge pour l'affectation.
*\param S Solide
*\return Solide
*/
Solide & Solide:: operator=(const Solide &S){
	
	assert(this != &S);
	solide.resize(S.solide.size());
	for(int i=0; i<S.solide.size(); i++){
		solide[i]= S.solide[i];
	}
}
/*!
*\fn void Solide::Affiche()
*\brief Fonction auxiliaire utile pour les tests.
*/
void Solide::Affiche(){
	
  for(int i=0; i<solide.size(); i++){
		cout<<"Particule "<<i<<endl;
    solide[i].Affiche();
  }

}
/*!
*\fn void Solide::Init(const char* s)
*\brief Initialisation du solide &agrave; partir d'un fichier. 
*\param s maillage solide
*\return void
*/
void Solide::Init(const char* s, const bool& rep, const int& numrep, const double& rho){
  std::ifstream maillage(s,ios::in);
  if(maillage){
    // cout <<"ouverture de xt.vtk reussie" << endl;
  } else {
    cout <<"ouverture de " << s << " ratee" << endl;
  }

  //Recuperation du maillage solide
  int Npoint;
  string sp;
  maillage >> sp >> Npoint;
  const int nb_points = Npoint;
  
  vector<Point_3> Points(nb_points);
  
  for(int i=0;i<nb_points;i++){
    double x,y,z;
    maillage >> x >> y >> z;
    Points[i] = Point_3(x,y,z);
  }
  
  int Npart;
  string sP;
  maillage >> sP >> Npart;
  const int nb_particule = Npart;
  
  vector<Particule> P(nb_particule);
  
  //bool points_particules[nb_points][nb_particule];
  
  // for(int i=0;i<nb_particule;i++){
  //   for(int j=0;j<nb_points;j++){
  //     //Remise a zero du tableau
  //     points_particules[j][i] = false;
  //   }
  // }
  vector<int> particules_vertex[nb_points];
  
  for(int i=0;i<nb_particule;i++){
    int Nfaces;
    int fixe;
    double X,Y,Z,u,v,w,theta,phi,psi,xmin,ymin,zmin,xmax,ymax,zmax;
    string s;
    maillage >> s >> Nfaces >> fixe;
    maillage >> s >> X >> Y >> Z;
    Point_3 centre(X,Y,Z);
    maillage >> s >> u >> v >> w;
    maillage >> s >> theta >> phi >> psi;
    xmin = X;
    ymin = Y;
    zmin = Z;
    xmax = X;
    ymax = Y;
    zmax = Z;
    const int nb_faces = Nfaces;
    std::vector<Face> Faces(nb_faces);
    std::vector<Point_3> points_c;
    for(int j=0;j<nb_faces;j++){
      int Nvertex;
      maillage >> Nvertex;
      const int nb_vertex = Nvertex;
      std::vector<Vertex> Vertex(nb_vertex);
      for(int k=0;k<nb_vertex;k++){
				int p;
				maillage >> p;
				Vertex[k].pos = Points[p];
				points_c.push_back(Points[p]);
				Vertex[k].num = p;
				//points_particules[p][i] = true;
				particules_vertex[p].push_back(i);
				double x = (Points[p][0]);
				double y = (Points[p][1]);
				double z = (Points[p][2]);
				xmin = min(x,xmin);
				xmax = max(x,xmax);
				ymin = min(y,ymin);
				ymax = max(y,ymax);
				zmin = min(z,zmin);
				zmax = max(z,zmax);
      }
      int voisin;
      maillage >> voisin;
      Faces[j] = Face(Vertex, voisin);
    }
    
		Point_3 center_part= centroid(points_c.begin(), points_c.end());
	
		if(fixe==0 || fixe==1){
		  P[i] = Particule(center_part, xmin, ymin, zmin, xmax, ymax, zmax, Faces);
		} else {
		  P[i] = Particule(centre, xmin, ymin, zmin, xmax, ymax, zmax, Faces);
		}
		P[i].fixe = fixe;
    P[i].u = Vector_3(u,v,w);
    P[i].omega = Vector_3(theta,phi,psi);
		P[i].u_half = Vector_3(u,v,w);
		P[i].omega_half = Vector_3(theta,phi,psi);
  }
  //Boucle de mise a jour des particules sur les sommets du maillage
  //Mise a jour des distances a l'equilibre entre particules en meme temps
  for(int i=0;i<P.size();i++){
    for(int j=0;j<P[i].faces.size();j++){
      for(int k=0;k<P[i].faces[j].size();k++){
	for(vector<int>::iterator l=particules_vertex[P[i].faces[j].vertex[k].num].begin();l!=particules_vertex[P[i].faces[j].vertex[k].num].end();l++){
	  //if(points_particules[P[i].faces[j].vertex[k].num][l]){
	    P[i].faces[j].vertex[k].particules.push_back(*l);
	    //cout << i << " " << j << " " << k << " " <<  P[i].faces[j].vertex[k].num << " " << l << endl;
			//cout << P[i].faces[j].vertex[k].num << " " << l << endl;
	    //getchar();
	    //}
	}
	if(P[i].faces[j].voisin>=0){
	  P[i].faces[j].D0 = sqrt((squared_distance(P[i].x0,P[P[i].faces[j].voisin].x0)));
	}
      }
    }
  }

  for(int i=0; i<P.size(); i++){
    solide.push_back(P[i]);
  }
  
  //Initialisation de la position et de l'inertie du solide
  for(int i=0; i<solide.size(); i++){
    solide[i].Dx = Vector_3(0.,0.,0.);
    solide[i].Dxprev = Vector_3(0.,0.,0.);
    solide[i].Fi = Vector_3(0.,0.,0.);
    solide[i].Ff = Vector_3(0.,0.,0.);
    solide[i].Ffprev = Vector_3(0.,0.,0.);
    solide[i].Mi = Vector_3(0.,0.,0.);
    solide[i].Mf = Vector_3(0.,0.,0.);
    solide[i].Mfprev = Vector_3(0.,0.,0.);
    solide[i].e = Vector_3(0.,0.,0.);
    solide[i].eprev = Vector_3(0.,0.,0.);
    solide[i].Inertie(rho);
    solide[i].mvt_t = Aff_transformation_3(1,0,0,0,1,0,0,0,1);
    solide[i].mvt_tprev = Aff_transformation_3(1,0,0,0,1,0,0,0,1);
  }

  //En cas de reprise
  if(rep){
    std::ostringstream oss;
    oss << "resultats/solide" << numrep << ".vtk";
    string s = oss.str();
    const char* nom = s.c_str();
    std::ifstream init(nom,std::ios::in);
    string dump;
    int nb_part = solide.size();
    int nb_triangles = 0.;
    for(int it=0; it<nb_part; it++){
      nb_triangles += solide[it].triangles.size();
    }
    for(int it=0;it<5*nb_triangles+13;it++){
      getline(init,dump);
    }
    cout << dump << endl;
    //Recuperation du deplacement
    for(int it=0; it<nb_part; it++){
      double Dx,Dy,Dz;
      init >> Dx >> Dy >> Dz;
      solide[it].Dx = Vector_3(Dx,Dy,Dz);
      for(int l= 0; l<solide[it].triangles.size(); l++){
	getline(init,dump);
      }
    }
    getline(init,dump);
    getline(init,dump);
    cout << dump << endl;
    //Recuperation de la vitesse
    for(int it=0; it<nb_part; it++){
      double u,v,w;
      init >> u >> v >> w;
      solide[it].u = Vector_3(u,v,w);
      for(int l= 0; l<solide[it].triangles.size(); l++){
	getline(init,dump);
      }
    }
    getline(init,dump);
    getline(init,dump);
    cout << dump << endl;
    //Recuperation du vecteur de rotation
    for(int it=0; it<nb_part; it++){
      double ex,ey,ez;
      init >> ex >> ey >> ez;
      solide[it].e = Vector_3(ex,ey,ez);
      for(int l= 0; l<solide[it].triangles.size(); l++){
	getline(init,dump);
      }
    }
    getline(init,dump);
    getline(init,dump);
    cout << dump << endl;
    //Recuperation de la vitesse de rotation
    for(int it=0; it<nb_part; it++){
      double omegax,omegay,omegaz;
      init >> omegax >> omegay >> omegaz;
      solide[it].omega = Vector_3(omegax,omegay,omegaz);
      for(int l= 0; l<solide[it].triangles.size(); l++){
	getline(init,dump);
      }
    }
    init.close();
    //Mise a jour de differents parametres
    for(int i=0; i<solide.size(); i++){
      solide[i].Dxprev = solide[i].Dx;
      solide[i].eprev = solide[i].e;
      double rot[3][3];
      //Recuperation de la matrice de rotation
      double e0 = sqrt(abs(1.-(solide[i].e.squared_length())));
      rot[0][0] = 1.-2.*(solide[i].e[1]*solide[i].e[1]+solide[i].e[2]*solide[i].e[2]);
      rot[0][1] = 2.*(-e0*solide[i].e[2]+solide[i].e[0]*solide[i].e[1]);
      rot[0][2] = 2.*(e0*solide[i].e[1]+solide[i].e[0]*solide[i].e[2]);
      rot[1][0] = 2.*(e0*solide[i].e[2]+solide[i].e[1]*solide[i].e[0]);
      rot[1][1] = 1.-2.*(solide[i].e[0]*solide[i].e[0]+solide[i].e[2]*solide[i].e[2]);
      rot[1][2] = 2.*(-e0*solide[i].e[0]+solide[i].e[1]*solide[i].e[2]);
      rot[2][0] = 2.*(-e0*solide[i].e[1]+solide[i].e[2]*solide[i].e[0]);
      rot[2][1] = 2.*(e0*solide[i].e[0]+solide[i].e[2]*solide[i].e[1]);
      rot[2][2] = 1.-2.*(solide[i].e[0]*solide[i].e[0]+solide[i].e[1]*solide[i].e[1]);
      Aff_transformation_3 rotation(rot[0][0],rot[0][1],rot[0][2],rot[1][0],rot[1][1],rot[1][2],rot[2][0],rot[2][1],rot[2][2]);
      Aff_transformation_3 translation(Vector_3(Point_3(0.,0.,0.),solide[i].x0)+solide[i].Dx);

      Aff_transformation_3 translation_inv(Vector_3(solide[i].x0,Point_3(0.,0.,0.)));
      solide[i].mvt_tprev = solide[i].mvt_t;
      solide[i].mvt_t = translation*(rotation*translation_inv);
    }
    update_triangles();
  }
  
}

/*!
*\fn void Solide::Solve_position(double dt, bool flag_2d)
*\brief Mise &agrave; jour de la position du solide.
*\param dt pas de temps
*\warning <b> Proc&eacute;dure sp&eacute;cifique au solide! </b>
*\return void
*/
void Solide::Solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  for(int i=0;i<size();i++){
    solide[i].solve_position(dt, flag_2d, t, T);
  }
  //breaking_criterion();
  update_triangles();
	for(int i=0;i<size();i++){
	  //double x_min = solide[i].max_x, y_min=solide[i].max_y, z_min=solide[i].max_z, x_max =solide[i].max_x, y_max=solide[i].max_y, z_max =solide[i].max_z;
	
	  for(std::vector<Triangle_3>::iterator it=solide[i].triangles.begin();it!=solide[i].triangles.end();it++){
	    for(int k=0;k<3;k++){
	      solide[i].bbox = Bbox(min(solide[i].bbox.xmin(),((*it).vertex(k).x())),min(solide[i].bbox.ymin(),((*it).vertex(k).y())),min(solide[i].bbox.zmin(),((*it).vertex(k).z())),max(solide[i].bbox.xmax(),((*it).vertex(k).x())),max(solide[i].bbox.ymax(),((*it).vertex(k).y())),max(solide[i].bbox.zmax(),((*it).vertex(k).z())));
	    }
	  }
	  
	  /*for(int j=0;j<solide[i].triangles.size();j++){
			
	    
			x_min = min( x_min ,min((solide[i].triangles[j][0].x()), min((solide[i].triangles[j][1].x()), (solide[i].triangles[j][2].x()) )) );
			
			y_min = min( y_min ,min((solide[i].triangles[j][0].y()), min((solide[i].triangles[j][1].y()), (solide[i].triangles[j][2].y()) )) );
			
			z_min = min( z_min ,min((solide[i].triangles[j][0].z()), min((solide[i].triangles[j][1].z()), (solide[i].triangles[j][2].z()) )) );
			
			x_max = max( x_max ,max((solide[i].triangles[j][0].x()), max((solide[i].triangles[j][1].x()), (solide[i].triangles[j][2].x()) )) );
			
			y_max = max( y_max ,max((solide[i].triangles[j][0].y()), max((solide[i].triangles[j][1].y()), (solide[i].triangles[j][2].y()) )) );
			
			z_max = max( z_max ,max((solide[i].triangles[j][0].z()), max((solide[i].triangles[j][1].z()), (solide[i].triangles[j][2].z()) )) );
		}
		solide[i].min_x = x_min; solide[i].min_y = y_min; solide[i].min_z = z_min;
		solide[i].max_x = x_max; solide[i].max_y = y_max; solide[i].max_z = z_max;*/
	}
	
}

/*void Solide::stock_def_plastique(const double &dt) {
  for(int i=0;i<size();i++){
    Particule* P = &solide[i];
    for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
      if((*F).voisin>=0 && (*F).plastifie){
	int part = (*F).voisin;
	double S = (*F).S;
	//Demi-incrément car on passe 2 fois par face...
	(*F).def_plas_cumulee += 0.5 * abs( ((*P).u - solide[part].u) * ((*P).u - solide[part].u) ) / (*F).D0 * dt; //Bien codé ? Semblerait pas...
      }
    }
  }
  }*/

/*!
*\fn void Solide::Solve_vitesse(double dt)
*\brief Calcul de la vitesse du solide.
*\param dt pas de temps
*\warning <b> Proc&eacute;dure sp&eacute;cifique au solide! </b>
*\return void
*/
void Solide::Solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  for(int i=0;i<size();i++){
    solide[i].solve_vitesse(dt, flag_2d, Amort, t , T);
  }
}

/*void Solide::Solve_vitesse_plas(const double& dt, const bool& flag_2d){
  for(int i=0;i<size();i++){
    solide[i].solve_vitesse_plas(dt, flag_2d);
  }
  }*/

/*!
*\fn void Solide::Forces(int N_dim, double nu, double E)
*\brief Calcul des forces. 
*\warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
*\return void
*/
void Solide::Forces(const int& N_dim, const double& nu, const double& E, const double& dt, const double& t, const double& T){
  Forces_internes(N_dim,nu,E, dt);
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (*P).Fi = (*P).Fi + Forces_externes((*P).x0+(*P).Dx,(*P).e,t,T);
    //(*P).Fi_plas = (*P).Fi_plas + Forces_externes((*P).x0+(*P).Dx,(*P).e);
    (*P).Mi = (*P).Mi + Moments_externes((*P).x0+(*P).Dx,(*P).e);
  }
}



/*!
*\fn void Solide::Forces_internes(int N_dim, double nu, double E)
*\brief Calcul des forces internes. 
*\warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
*\return void
*/

// double Face::Forces_elas(Particule* P, std::vector<Particule> solide, const double &nu, const double &E) {
//   //Mettre ici forces élastiques !
//   int part = voisin;
//   Point_3 xi((*P).x0); //Position particule i en config initiale
//   Point_3 xj(solide[part].x0); //Position particule j en config initiale
//   Vector_3 lij(xi, xj);
//   //double dij = sqrt( lij.squared_length() );
//   Vector_3 nIJ = lij / D0;
//   double Dij_n = ((*P).Dx - solide[part].Dx ) * nIJ;
//   //double K = E/(1.-2.*nu)/3.;
//   return S/6.*E*(Dij_n/D0 - def_plas_cumulee); //Force élastique du lien
// }

// double Face::Forces_plas(Particule* P, std::vector<Particule> solide, const double &n, const double &B) {
//   //Mettre ici forces élastiques !
//   int part = voisin;
//   Point_3 xi((*P).x0); //Position particule i en config initiale
//   Point_3 xj(solide[part].x0); //Position particule j en config initiale
//   Vector_3 lij(xi, xj);
//   //double dij = sqrt( lij.squared_length() );
//   Vector_3 nIJ = lij / D0;
//   double Dij_n = ((*P).Dx - solide[part].Dx /*- cross_product((*P).omega + solide[part].omega, lij / 2.)*/ ) * nIJ;
//   double volume_diam = D0 / 2. * S / 3.;
//   //return -signe(Dij_n)* B * volume_diam / pow(D0, n+1) * pow(abs(Dij_n), n) * nIJ; //Force plastique du lien
//   return -signe(Dij_n) * B * volume_diam * pow(def_plas_cumulee, n); //Force plastique du lien
// }

void Solide::Forces_internes(const int& N_dim, const double& nu, const double& E, const double& dt){

  bool plastifie = false;
  
  //Calcul de la contrainte dans chaque particule
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    //(*P).Volume_libre();
    (*P).discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    (*P).discrete_gradient.col2 = Vector_3(0., 0., 0.);
    (*P).discrete_gradient.col3 = Vector_3(0., 0., 0.);
    for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
      if((*F).voisin>=0){
	int part = (*F).voisin;
	/*Point_3 xi((*P).x0); //Position particule i en config initiale
        Point_3 xj(solide[part].x0); //Position particule j en config initiale
	Vector_3 lij(xi, xj);
	Vector_3 nIJ = lij / (*F).D0;*/
	Vector_3 nIJ = (*F).normale;
	//double Dij_n = (solide[part].Dx - (*P).Dx ) * nIJ; //Quadrature au point gauche
	Matrix Dij_n(tens_sym(solide[part].Dx + solide[part].u * dt/2. - (*P).Dx - (*P).u * dt/2.,  nIJ) ); //Quadrature au point milieu pour calcul des forces !
	(*P).discrete_gradient += (*F).S / 2. * Dij_n / (*P).volume();

      }
    }
    //cout << "Trace dev Def : " << (((*P).discrete_gradient).dev()).tr() << endl;
    (*P).contrainte = lambda * ((*P).discrete_gradient - (*P).epsilon_p).tr() * unit() + 2*mu * ((*P).discrete_gradient - (*P).epsilon_p);
    //cout << "Trace dev Contrainte : " << (((*P).contrainte).dev()).tr() << endl;

    //Mettre ces valeurs dans le param.dat !!!!!
    double B = 292000000.; //En Pa. JC.
    double n = .31; //JC.
    double A = 90000000.; //En Pa. Vient de JC
    double H = 0.; //60000000000.; //En Pa. Moitié du module de Young
	
    (*P).seuil_elas = A; // + B * pow((*P).def_plas_cumulee, n);

    if((P->contrainte - H * P->epsilon_p).VM() > (*P).seuil_elas) { //On sort du domaine élastique.
      plastifie = true;
      //Matrix n_elas(((*P).contrainte).dev() / (((*P).contrainte).dev()).norme() ); //Normale au domaine élastique de Von Mises
      Matrix n_elas( 1. / (((*P).contrainte).dev()).norme() * ((*P).contrainte).dev() ); //Normale au domaine élastique de Von Mises
      /*if((*P).n_elas_prev == -n_elas)
	cout << "Chargement dans sens oppose !" << endl;
      (*P).n_elas_prev = n_elas;*/
      //cout << "Trace n_elas : " << n_elas.tr() << endl;
      //cout << "Norme n_elas : " << n_elas.norme() << endl;
      //double delta_p = pow(((*P).contrainte.VM() - A) / B, 1./n) - (*P).def_plas_cumulee;
      double delta_p = ((P->contrainte - H * P->epsilon_p).VM() - A) / (2*mu + H);
      //(*P).def_plas_cumulee = pow(((*P).contrainte.VM() - A) / B, 1./n); //Nouvelle déformation plastique.
      //cout << "Def plastique cumulee : " << (*P).def_plas_cumulee << endl;
      (*P).epsilon_p += delta_p * n_elas;
      //cout << "Trace def plas : " << ((*P).epsilon_p).tr() << endl; //Pb ! Non-nulle !!!!
      //cout << "Norme def plas : " << ((*P).epsilon_p).norme() << endl;
      
      //((*P).contrainte.VM() - A) / E * n_elas;  //* signe( (*P).contrainte ); //Plasticité parfaite
      //(*P).contrainte = A * signe( (*P).contrainte );
    }
  }

  if(plastifie)
    cout << "Plastification dans ce pas de temps !" << endl;

  
  //Calcul des forces pour chaque particule
  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (*P).Fi = Vector_3(0.,0.,0.);
    for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
      if((*F).voisin>=0){
	int part = (*F).voisin;
	/*Point_3 xi((*P).x0); //Position particule i en config initiale
        Point_3 xj(solide[part].x0); //Position particule j en config initiale
	Vector_3 lij(xi, xj);
	Vector_3 nIJ = lij / (*F).D0;*/
	Vector_3 nIJ = (*F).normale;

	//Il faut passer les contraintes dans le repère local de chaque face (nIJ, sIJ, tIJ) ????
        Vector_3 Fij_elas( (*F).S / 2. * ( ((*P).contrainte + solide[part].contrainte).tr() / 2. * nIJ + ((*P).contrainte + solide[part].contrainte) / 2. * nIJ ) ); //Force du lien IJ !
	//cout << "Force : " << Fij_elas << endl;

	(*P).Fi = (*P).Fi + Fij_elas; // * nIJ; //Force sur particule
	
	  
	
	//Force de rappel elastique
	//(*P).Fi = (*P).Fi + S/(*F).D0*E/(1.+nu)*Delta_u;
	//Force de deformation volumique
	//(*P).Fi = (*P).Fi + S*E*nu/(1.+nu)/(1.-2.*nu)*epsilonIJ*(nIJ+Delta_u/DIJ-(Delta_u*nIJ)/DIJ*nIJ); //Changer cette expression aussi ? Surement...
	//(*P).Fi = (*P).Fi + S*E/(1.-2.*nu)*epsilonIJ*nIJ;
	//Moment des forces appliquees
	/*(*P).Mi = (*P).Mi + cross_product((*P).mvt_t(XC1),S/(*F).D0*E/(1.+nu)*Delta_u);
	  (*P).Mi = (*P).Mi + cross_product((*P).mvt_t(XC1),S*E*nu/(1.+nu)/(1.-2.*nu)*epsilonIJ*(nIJ+Delta_u/DIJ-(Delta_u*nIJ)/DIJ*nIJ));
	  //Moments de flexion/torsion
	  double kappa = 1.;
	  double alphan = (2.+2.*nu-kappa)*E/4./(1.+nu)/S*((*F).Is+(*F).It);
	  double alphas = E/4./(1.+nu)/S*((2.+2.*nu+kappa)*(*F).Is-(2.+2.*nu-kappa)*(*F).It);
	  double alphat = E/4./(1.+nu)/S*((2.+2.*nu+kappa)*(*F).It-(2.+2.*nu-kappa)*(*F).Is);
	  (*P).Mi = (*P).Mi + S/(*F).D0*(alphan*cross_product((*P).mvt_t((*F).normale),solide[part].mvt_t((*F).normale))+alphas*cross_product((*P).mvt_t((*F).s),solide[part].mvt_t((*F).s))+alphat*cross_product((*P).mvt_t((*F).t),solide[part].mvt_t((*F).t))); */
	// cout << alphan << " " << alphas << " " << alphat  << endl;
	// cout << S/(*F).D0*alphan << endl;
	// cout << "D=" << (*F).D0 << " I=" << (*P).I[0] << " " << (*P).I[1] << " " << (*P).I[2] << " S=" << S << endl;
	// cout << "DT=" << sqrt(min(min((*P).I[0],(*P).I[1]),(*P).I[2])*(*F).D0/S/alphan) << endl;
	// getchar();
      }
    }
  }
}

/*!
*\fn double Solide::Energie(int N_dim, double nu, double E)
*\brief Calcul d'&eacute;nergie. 
*\warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
*\return void
*/
double Solide::Energie(const int& N_dim, const double& nu, const double& E){
  return Energie_cinetique()+Energie_potentielle(N_dim, nu, E);
}
/*!
*\fn double Solide::Energie_cinetique()
*\brief Calcul d'&eacute;nergie cin&eacute;tique. 
*\warning  <b> Proc&eacute;dure sp&eacute;cifique au solide! </b> 
*\return void
*/
double Solide::Energie_cinetique(){
  double E = 0.;
  for(int i=0;i<size();i++){
    double u2 = (solide[i].u.squared_length());
    E += 1./2.*solide[i].m*u2;
    //Calcul de -1/2*tr(D j(Q^T omega)) = 1/2*(I1*Omega1^2+I2*Omega2^2+I3*Omega3^2)
    double Q[3][3];
    double rot[3][3];
    double e0 = sqrt(abs(1.-(solide[i].e.squared_length())));
    //Recuperation de la matrice de rotation
    rot[0][0] = 1.-2.*(solide[i].e[1]*solide[i].e[1]+solide[i].e[2]*solide[i].e[2]);
    rot[0][1] = 2.*(-e0*solide[i].e[2]+solide[i].e[0]*solide[i].e[1]);
    rot[0][2] = 2.*(e0*solide[i].e[1]+solide[i].e[0]*solide[i].e[2]);
    rot[1][0] = 2.*(e0*solide[i].e[2]+solide[i].e[1]*solide[i].e[0]);
    rot[1][1] = 1.-2.*(solide[i].e[0]*solide[i].e[0]+solide[i].e[2]*solide[i].e[2]);
    rot[1][2] = 2.*(-e0*solide[i].e[0]+solide[i].e[1]*solide[i].e[2]);
    rot[2][0] = 2.*(-e0*solide[i].e[1]+solide[i].e[2]*solide[i].e[0]);
    rot[2][1] = 2.*(e0*solide[i].e[0]+solide[i].e[2]*solide[i].e[1]);
    rot[2][2] = 1.-2.*(solide[i].e[0]*solide[i].e[0]+solide[i].e[1]*solide[i].e[1]);
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	Q[j][k] = rot[j][0]*solide[i].rotref[0][k];
	Q[j][k] += rot[j][1]*solide[i].rotref[1][k];
	Q[j][k] += rot[j][2]*solide[i].rotref[2][k];
      }
    }  
    double Omega[3];
    Omega[0] = Omega[1] = Omega[2] = 0.;
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
	Omega[j] += (solide[i].omega[k]*Q[k][j]);
      }
    }
    E += 1./2.*(solide[i].I[0]*Omega[0]*Omega[0]+solide[i].I[1]*Omega[1]*Omega[1]+solide[i].I[2]*Omega[2]*Omega[2]);
  }
  return E;
}

double Solide::Energie_potentielle(const int& N_dim, const double& nu, const double& E){
  double Ep = 0.;

  double B = 292000000.; //En Pa. JC.
  double n = .31; //JC.
  double A = 90000000.; //En Pa. Vient de JC

  for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
    Ep += ( 0.5 * contraction_double((*P).contrainte, (*P).discrete_gradient - (*P).epsilon_p)  + B * pow((*P).def_plas_cumulee, 1. + n) / (n + 1.) + A * (*P).def_plas_cumulee ) * (*P).volume();
  }
  return Ep;
}

double Solide::pas_temps(const double& t, const double& T, const double& cfls, const double& E, const double& nu, const double& rhos){
  double eps = 1e-14;//std::numeric_limits<double>::epsilon();
  double dt = std::numeric_limits<double>::infinity();
  //Restriction CFL sur la vitesse de rotation
  for(int i=0;i<size();i++){
    double dt1 = cfls*0.26/(abs((solide[i].omega[0]))+abs((solide[i].omega[1]))+abs((solide[i].omega[2]))+eps);
    dt = min(dt,dt1); 
  }
  //Restriction CFL liee aux forces internes
  double cs = sqrt(E*(1.-nu)/rhos/(1.+nu)/(1.-2.*nu));
  //Calcul du rayon de la sphï¿½re inscrite
  double sigma = 100000.;
  for(int i=0;i<size();i++){
    for(int j=0;j<solide[i].faces.size();j++){
      sigma = min(sigma,solide[i].faces[j].D0);
    }
  }
  for(int i=0;i<size();i++){
    for(int j=0;j<solide[i].faces.size();j++){
      if(solide[i].faces[j].voisin>=0){
	dt = min(dt,cfls*solide[i].faces[j].D0/cs);
	//dt = min(dt,cfls*0.26*sqrt(pow(sigma,5)/solide[i].faces[j].S/solide[i].faces[j].D0)/cs);
	double Imin = min(min(solide[i].I[0],solide[i].I[1]),solide[i].I[2]);
	double S = solide[i].faces[j].S;
	double D0 = solide[i].faces[j].D0;
	double kappa = 1.;
	double alphan = (2.+2.*nu-kappa)*E/4./(1.+nu)/S*(solide[i].faces[j].Is+solide[i].faces[j].It);
	double alphas = E/4./(1.+nu)/S*((2.+2.*nu+kappa)*solide[i].faces[j].Is-(2.+2.*nu-kappa)*solide[i].faces[j].It);
	double alphat = E/4./(1.+nu)/S*((2.+2.*nu+kappa)*solide[i].faces[j].It-(2.+2.*nu-kappa)*solide[i].faces[j].Is);
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphan));
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphas));
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphat));
      }
    }
  }
  dt = min(dt,T-t);
  return dt;
}

/*!
*\fn void Solide::update_triangles()
*\brief Mise &agrave; jour de l'interface fluide - solide.
*\details Mise &agrave; jour des \a Particule.triangles_prev, \a Particule.triangles, \a Particule.normales_prev, \a Particule.normales, \a Particule.fluide_prev, \a Particule.fluide, \a Particule.Points_interface_prev, \a Particule.Points_interface, \a Particule.Triangles_interface_prev, \a Particule.Triangles_interface, \a Particule.Position_Triangles_interface_prev et \a Particule.Position_Triangles_interface.
*\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
*\return void
*/
void Solide::update_triangles(){
	for(int i=0;i<solide.size();i++){
		solide[i].triangles_prev = solide[i].triangles;
		solide[i].normales_prev = solide[i].normales;
		solide[i].fluide_prev = solide[i].fluide;
		for(int it=0;it<solide[i].triangles.size();it++){
			solide[i].Points_interface_prev[it] = solide[i].Points_interface[it];
			solide[i].Triangles_interface_prev[it] = solide[i].Triangles_interface[it];
			solide[i].Position_Triangles_interface_prev[it] = solide[i].Position_Triangles_interface[it];
			solide[i].Points_interface[it].erase(solide[i].Points_interface[it].begin(),solide[i].Points_interface[it].end());
			solide[i].Triangles_interface[it].erase(solide[i].Triangles_interface[it].begin(),solide[i].Triangles_interface[it].end());	solide[i].Position_Triangles_interface[it].erase(solide[i].Position_Triangles_interface[it].begin(),
                                                       solide[i].Position_Triangles_interface[it].end());
		}
		solide[i].triangles.erase(solide[i].triangles.begin(),solide[i].triangles.end());
		solide[i].normales.erase(solide[i].normales.begin(),solide[i].normales.end());
		solide[i].fluide.erase(solide[i].fluide.begin(),solide[i].fluide.end());
		solide[i].vide.erase(solide[i].vide.begin(),solide[i].vide.end());
		
		//Calcul de la nouvelle position des triangles
		for(int f=0;f<solide[i].faces.size();f++){
			Point_3 s,r,v,t;
			
			if(solide[i].faces[f].size() == 3){
				vector<Point_3> ri,vi,si ;
				for(int part=0; part<solide[i].faces[f].vertex[0].size();part++){
					int p = solide[i].faces[f].vertex[0].particules[part];
					ri.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[0].pos));
				}
				r = centroid(ri.begin(),ri.end());
				
				
				for(int part=0;part<solide[i].faces[f].vertex[1].size();part++){
					int p = solide[i].faces[f].vertex[1].particules[part];
					vi.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[1].pos));
				}
				v = centroid(vi.begin(),vi.end());

				for(int part=0;part<solide[i].faces[f].vertex[2].size();part++){
					int p = solide[i].faces[f].vertex[2].particules[part];
					si.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[2].pos));
				}
				s = centroid(si.begin(),si.end());
				
				Vector_3 vect0(r,v);
				Vector_3 vect1(r,s);
				Triangle_3 Tri(r,v,s);
				solide[i].triangles.push_back(Tri);
				Vector_3 normale = cross_product(vect0,vect1);
				normale = normale*(1./sqrt((normale.squared_length())));
				solide[i].normales.push_back(normale);
				if(solide[i].faces[f].voisin < 0){
					solide[i].fluide.push_back(true);
				} else {
					solide[i].fluide.push_back(false);
				}
				if( solide[i].faces[f].voisin == -2){
					solide[i].vide.push_back(true);
				} else {
					solide[i].vide.push_back(false);
				}
			}
// 	 else if(flag_2d){
//  cout<<"update tag 0"<<endl;
// 			 vector<Point_3> si,ri,vi,ti;
// 			 for(int part=0;part<solide[i].faces[f].vertex[0].size();part++){
// 				 int p = solide[i].faces[f].vertex[0].particules[part];
// 				 si.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[0].pos));
// 			 }
// 			 s = centroid(ri.begin(),ri.end());
// 			 for(int part=0;part<solide[i].faces[f].vertex[1].size();part++){
// 				 int p = solide[i].faces[f].vertex[1].particules[part];
// 				 ri.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[1].pos));
// 			 }
// 			 r = centroid(ri.begin(),ri.end());
// 			 for(int part=0;part<solide[i].faces[f].vertex[2].size();part++){
// 				 int p = solide[i].faces[f].vertex[2].particules[part];
// 				 vi.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[2].pos));
// 			 }
// 			 v = centroid(vi.begin(),vi.end());
// 			 
// 			 for(int part=0;part<solide[i].faces[f].vertex[3].size();part++){
// 				 int p = solide[i].faces[f].vertex[3].particules[part];
// 				 ti.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[3].pos));
// 			 }
// 			 t = centroid(vi.begin(),vi.end());
// 			 
// 			 
// 			 Vector_3 vect0(r,s);
// 			 Vector_3 vect1(r,v);
// 			 Triangle_3 Tri1(s,r,v);
// 			 solide[i].triangles.push_back(Tri1);
// 			 Vector_3 normale = CGAL::cross_product(vect0,vect1);
// 			 normale = normale*(1./sqrt((normale.squared_length())));
// 			 solide[i].normales.push_back(normale);
// 			 
// 			 Vector_3 vect2(v,s);
// 			 Vector_3 vect3(v,t);
// 			 Triangle_3 Tri2(v,t,s);
// 			 solide[i].triangles.push_back(Tri2);
// 			 Vector_3 normale2 = CGAL::cross_product(vect2,vect3);
// 			 normale2 = normale2*(1./sqrt((normale2.squared_length())));
// 			 solide[i].normales2.push_back(normale2);
// 			 
// 			 
// 			 
// 			 
// 			 if(solide[i].faces[f].voisin < 0){
// 				 solide[i].fluide.push_back(true);
// 				 solide[i].fluide.push_back(true);
// 			 } 
// 			 else {
// 				 solide[i].fluide.push_back(false);
// 				 solide[i].fluide.push_back(false);
// 			 }
// 			 if( solide[i].faces[f].voisin == -2){
// 				 solide[i].vide.push_back(true);
// 				 solide[i].vide.push_back(true);
// 			 } else {
// 				 solide[i].vide.push_back(false);
// 				 solide[i].vide.push_back(false);
// 			 }
// 			 cout<<"update tag 1"<<endl;
// 	 }
		else{
			
			vector<Point_3> si;
			si.push_back(solide[i].mvt_t(solide[i].faces[f].centre));
			int j = solide[i].faces[f].voisin;
			if(j>=0){
				si.push_back(solide[j].mvt_t(solide[i].faces[f].centre));
			}
			s = centroid(si.begin(),si.end());
			
			for(int k=0;k<solide[i].faces[f].size();k++){
			  int kp = (k+1)%(solide[i].faces[f].size());
			  vector<Point_3> ri,vi;
			  for(int part=0;part<solide[i].faces[f].vertex[k].size();part++){
			    int p = solide[i].faces[f].vertex[k].particules[part];
			    ri.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[k].pos));
			  }
			  r = centroid(ri.begin(),ri.end());
			  for(int part=0;part<solide[i].faces[f].vertex[kp].size();part++){
			    int p = solide[i].faces[f].vertex[kp].particules[part];
			    vi.push_back(solide[p].mvt_t(solide[i].faces[f].vertex[kp].pos));
			  }
			  v = centroid(vi.begin(),vi.end());
			  Vector_3 vect0(s,r);
			  Vector_3 vect1(s,v);
			  Triangle_3 Tri(s,r,v);
			  solide[i].triangles.push_back(Tri);
			  Vector_3 normale = cross_product(vect0,vect1);
			  normale = normale*(1./sqrt((normale.squared_length())));			  
			  solide[i].normales.push_back(normale);
			  if(solide[i].faces[f].voisin < 0){
			    solide[i].fluide.push_back(true);
			  } 
			  else {
			    solide[i].fluide.push_back(false);
			  }
			  if( solide[i].faces[f].voisin == -2){
			    solide[i].vide.push_back(true);
			  } else {
			    solide[i].vide.push_back(false);
			  }
			}
		}
			
		}//Calcul de la nouvelle position des triangles
		
	}
}


// /*!
// *\fn double Error(Solide& S1, Solide& S2)
// *\brief Calcul d'erreur.
// *\details Fonction appell&eacute;e dans le cas d'un sch&eacute;ma semi-implicite.  Crit&egrave;re d'arr&ecirc;t:  \n
// \f{eqnarray*}{ error = max( \, \Vert S1.solide[i].Dx -  S2.solide[i].Dx \, \Vert_{\infty} + h_{max} \Vert \, S1.solide[i].e -  S2.solide[i].e \, \Vert_{\infty})_i  \f}\n
// \f{eqnarray*}{
// 	h_{max}=& max(abs(S1.max\_x - S1.min\_x),abs(S1.max\_y - S1.min\_y), \\ 
// 	& abs(S1.max\_z - S1.min\_z), abs(S2.max\_x - S2.min\_x),\\ 
// 	& abs(S2.max\_y - S2.min\_y),abs(S2.max\_z - S2.min\_z)).
// 				 \f}
				 
// *\param S1 \a Solide au temps t+k
// *\param S2 \a Solide au temps t+k-1
// *\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
// *\return double
// */

// double Error(Solide& S1, Solide& S2){
	
// 	double erreur = -1.;
	
// 	for(int it=0; it<S1.size(); it++){
		
// 	  double h_max1 = std::max(std::max((S1.solide[it].bbox.xmax() - S1.solide[it].bbox.xmin()),(S1.solide[it].bbox.ymax() - S1.solide[it].bbox.ymin())),              (S1.solide[it].bbox.zmax() - S1.solide[it].bbox.zmin())); 
// 	  double h_max2 = std::max(std::max((S2.solide[it].bbox.xmax() - S2.solide[it].bbox.xmin()),(S2.solide[it].bbox.ymax() - S2.solide[it].bbox.ymin())),              (S2.solide[it].bbox.zmax() - S2.solide[it].bbox.zmin())); 
// 	  double h_max = max(h_max1, h_max2);
// 	  double err1 = std::max(std::max(abs((S1.solide[it].Dx[0] - S2.solide[it].Dx[0])), abs((S1.solide[it].Dx[1] - S2.solide[it].Dx[1]) )), abs((S1.solide[it].Dx[2] - S2.solide[it].Dx[2]))); 
// 	  double err2 = std::max(std::max(abs((S1.solide[it].e[0] - S2.solide[it].e[0])), abs((S1.solide[it].e[1] - S2.solide[it].e[1]))), abs((S1.solide[it].e[2] - S2.solide[it].e[2]))); ;
// 	  double erreur_temp = err1 + h_max * err2;
// 	  erreur = std::max(erreur_temp, erreur);
// 	}
	
// 	return erreur;
// }	

// /*!
// * \fn void Copy_f_m(Solide& S1, Solide& S2)
// *  \brief On copie les valeurs \a Ff et \a Mf du S2 dans S1.
// *  \details Fonction appell&eacute;e dans le cas d'un sch&eacute;ma semi-implicite. 
// *	\param S1 \a Solide au temps t
// *	\param S2 \a Solide au temps t+k
// *	\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
// *	\return void
// 	*/
// void Copy_f_m(Solide& S1, Solide& S2){
	
// 	for(int it=0; it<S1.size(); it++){
// 		S1.solide[it].Ff =  S2.solide[it].Ff ;
// 		S1.solide[it].Mf =  S2.solide[it].Mf ;
// 	}
	
// }	
// /*!
// * \fn bool inside_box(const Bbox& cell, const Point_3& P)
// *\brief Fonction qui renvoie true si P est dans Box et false sinon.
// *\param cell \a Box 
// *\param P \a un point
// *\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
// *\return bool
// */
// bool inside_box(const Bbox& cell, const Point_3& P){
	
//   /*bool in = false;
	
// 	if((cell.xmin() - P.x())<= eps_relat && (cell.ymin() - P.y())<= eps_relat &&
// 		(cell.zmin() - P.z())<= eps_relat && (cell.xmax() - P.x())>=-eps_relat &&
// 		(cell.ymax() - P.y())>=-eps_relat && (cell.zmax() - P.z())>=-eps_relat )
// 	{ in = true; }
	
// 	return in;*/
//   return ((cell.xmin() - P.x())<= eps_relat && (cell.ymin() - P.y())<= eps_relat && (cell.zmin() - P.z())<= eps_relat && (cell.xmax() - P.x())>=-eps_relat && (cell.ymax() - P.y())>=-eps_relat && (cell.zmax() - P.z())>=-eps_relat );
  
// }	

// /*!
// * \fn bool inside_convex_polygon(const Particule& S, const Point_3& P)
// *\brief Fonction qui renvoie true si P(point) est dans S(polygon convex) et false sinon.
// *\param S \a Particule
// *\param P \a un point
// *\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
// *\return bool
// */
// bool inside_convex_polygon(const Particule& S, const Point_3& P){
	
// 	bool in = false;
	
// 	//if((S.min_x - P.x())<= eps_relat && (S.min_y - P.y())<= eps_relat && (S.min_z - P.z())<= eps_relat && (S.max_x - P.x())>=-eps_relat && (S.max_y - P.y())>=-eps_relat && (S.max_z - P.z())>=-eps_relat )
// 	if(CGAL::do_overlap(S.bbox,P.bbox()))
// 	{
// 		if(S.cube) {in = true;}
// 		else{
// 			in = true;
// 			for(int l= 0; l<S.triangles.size() && in; l++){
// 			  const Point_3& vertex = S.triangles[l][0];
// 			  Vector_3 vect(P,vertex);
// 			  if(((vect*S.normales[l])) < 0.){in = false;}
// 			}
// 		}
// 	}
	
// 	return in;
// }	



// /*!
// *\fn bool box_inside_convex_polygon(const Particule& S, const Bbox& cell)
// *\brief Fonction qui renvoie true si cell(Box) est dans S(polygon convex) et false sinon.
// *\param S \a Particule
// *\param cell \a un Box
// *\warning <b> Proc&eacute;dure sp&eacute;cifique au couplage! </b>
// *\return bool
// */
// bool box_inside_convex_polygon(const Particule& S, const Bbox& cell){
	
// 	bool in = false;
	
// 	//if ((S.min_x - cell.xmin()) <= eps_relat && (S.min_y - cell.ymin() <= eps_relat) && (S.min_z - cell.zmin()) <= eps_relat && (S.max_x - cell.xmax() >=-eps_relat) && (S.max_y - cell.ymax()) >=-eps_relat && (S.max_z - cell.zmax()>=-eps_relat) ) 
// 	if(CGAL::do_overlap(S.bbox,cell))
// 	{
// 		if(S.cube) { return S.cube;}
		
// 		else{
			
// 			in = true;
			
// 			Point_3 p1(cell.xmin(),cell.ymin(),cell.zmin());
// 			in = inside_convex_polygon(S,p1);
// 			if(!in) {return in;}
			
// 			Point_3 p2(cell.xmax(),cell.ymax(),cell.zmax());
// 			in = inside_convex_polygon(S,p2);
// 			if(!in) {return in;}
			
// 			Point_3 p3(cell.xmax(),cell.ymin(),cell.zmin());
// 			in = inside_convex_polygon(S,p3);
// 			if(!in) {return in;}
			
// 			Point_3 p4(cell.xmax(),cell.ymin(),cell.zmax());
// 			in = inside_convex_polygon(S,p4);
// 			if(!in) {return in;}
			
// 			Point_3 p5(cell.xmax(),cell.ymax(),cell.zmin());
// 			in = inside_convex_polygon(S,p5);
// 			if(!in) {return in;}
			
// 			Point_3 p6(cell.xmin(),cell.ymax(),cell.zmin());
// 			in = inside_convex_polygon(S,p6);
// 			if(!in) {return in;}
			
// 			Point_3 p7(cell.xmin(),cell.ymax(),cell.zmax());
// 			in = inside_convex_polygon(S,p7);
// 			if(!in) {return in;}
			
// 			Point_3 p8(cell.xmin(),cell.ymin(),cell.zmax());
// 			in = inside_convex_polygon(S,p8);
// 			if(!in) {return in;}
			
// 		}
// 	}
	
// 	return in;
// }


// bool box_inside_tetra(const Tetrahedron &tetra, const Bbox& cell){
	
// 	bool in = false;
	
// 	//Bbox box_tetra= tetra.bbox();
	
// 	//if ((box_tetra.xmin() - cell.xmin()) <= eps_relat && (box_tetra.ymin() - cell.ymin() <= eps_relat) && (box_tetra.zmin() - cell.zmin()) <= eps_relat && (box_tetra.xmax() - cell.xmax() >=-eps_relat) && (box_tetra.ymax() - cell.ymax()) >=-eps_relat && (box_tetra.zmax() - cell.zmax()>=-eps_relat) )
// 	if(CGAL::do_overlap(tetra.bbox(),cell))
// 	{
// 			in = true;
			
// 			Point_3 p1(cell.xmin(),cell.ymin(),cell.zmin());
// 			in = inside_tetra(tetra,p1);
// 			if(!in) {return in;}
			
// 			Point_3 p2(cell.xmax(),cell.ymax(),cell.zmax());
// 			in = inside_tetra(tetra,p2);
// 			if(!in) {return in;}
			
// 			Point_3 p3(cell.xmax(),cell.ymin(),cell.zmin());
// 			in = inside_tetra(tetra,p3);
// 			if(!in) {return in;}
			
// 			Point_3 p4(cell.xmax(),cell.ymin(),cell.zmax());
// 			in = inside_tetra(tetra,p4);
// 			if(!in) {return in;}
			
// 			Point_3 p5(cell.xmax(),cell.ymax(),cell.zmin());
// 			in = inside_tetra(tetra,p5);
// 			if(!in) {return in;}
			
// 			Point_3 p6(cell.xmin(),cell.ymax(),cell.zmin());
// 			in = inside_tetra(tetra,p6);
// 			if(!in) {return in;}
			
// 			Point_3 p7(cell.xmin(),cell.ymax(),cell.zmax());
// 			in = inside_tetra(tetra,p7);
// 			if(!in) {return in;}
			
// 			Point_3 p8(cell.xmin(),cell.ymin(),cell.zmax());
// 			in = inside_tetra(tetra,p8);
// 			if(!in) {return in;}

// 	}
	
// 	return in;
// }



/*!
*\fn void Solide::Impression(int n)
*\brief Impression des r&eacute;sultats. 
*\param n num&eacute;ro de l'iteration en temps
*\return void
*/
void Solide::Impression(const int &n, const bool &reconstruction){ //Sortie au format vtk
  int nb_part = solide.size();
//Version avec reconstruction
  if(reconstruction){
    int nb_triangles = 0.;
    for(int it=0; it<nb_part; it++){
      nb_triangles += solide[it].triangles.size();
    }

//const char* solidevtk;
//{
    std::ostringstream oss;
    oss << "resultats/solide" << n << ".vtk";
    string s = oss.str();
    //cout << s << endl;
    const char* const solidevtk = s.c_str();
    //}
	
    //Ouverture des flux en donne en ecriture
    std::ofstream vtk;
    vtk.open(solidevtk,ios::out);
    if(vtk.is_open())
    {
      // cout <<"ouverture de xt.vtk reussie" << endl;
    } else {
      cout <<"ouverture de solide" << n << ".vtk rate" << endl;
    }
    vtk << setprecision(15);
    //Initialisation du fichier vtk
    vtk << "# vtk DataFile Version 3.0" << endl;
    vtk << "#Simulation Euler" << endl;
    vtk << "ASCII" << endl;
    vtk<<"\n";
    vtk << "DATASET UNSTRUCTURED_GRID" << endl;
    vtk << "POINTS " << 3*nb_triangles << " DOUBLE" << endl;
	
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++){
	vtk << solide[it].triangles[l][0][0] << " " << solide[it].triangles[l][0][1] << " " << solide[it].triangles[l][0][2] << endl;
	vtk << solide[it].triangles[l][1][0] << " " << solide[it].triangles[l][1][1] << " " << solide[it].triangles[l][1][2] << endl;
	vtk << solide[it].triangles[l][2][0] << " " << solide[it].triangles[l][2][1] << " " << solide[it].triangles[l][2][2] << endl;
      }
    }
    vtk<<"\n";
    vtk << "CELLS " << nb_triangles << " " << 4*nb_triangles<< endl;
    int num=0;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++){
	vtk << 3 << " " << 3*num << " " << 3*num+1 << " " << 3*num+2 << endl;
	num++;
      }
    }
    vtk << "\n";
    vtk << "CELL_TYPES " << nb_triangles << endl;
    for(int l= 0; l<nb_triangles; l++)
    {
      vtk << 5 << endl;
    }
    vtk << "\n";
    vtk << "CELL_DATA " << nb_triangles << endl;
    //Deplacement
    vtk << "VECTORS deplacement double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++)
      {
	vtk << solide[it].Dx[0] << " " << solide[it].Dx[1] << " " << solide[it].Dx[2] << endl;
      }
    }
    vtk << "\n";
    //Vitesse
    vtk << "VECTORS vitesse double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++)
      {
	vtk << solide[it].u[0] << " " << solide[it].u[1] << " " << solide[it].u[2] << endl;
      }
    }
    vtk << "\n";
    //Contrainte
    vtk << "TENSORS contraintes double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++)
      {
	vtk << solide[it].contrainte.col1[0] << " " << solide[it].contrainte.col1[1] << " " << solide[it].contrainte.col1[2] << endl;
	vtk << solide[it].contrainte.col2[0] << " " << solide[it].contrainte.col2[1] << " " << solide[it].contrainte.col2[2] << endl;
	vtk << solide[it].contrainte.col3[0] << " " << solide[it].contrainte.col3[1] << " " << solide[it].contrainte.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Déformations
    vtk << "TENSORS deformations double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++)
      {
	vtk << solide[it].discrete_gradient.col1[0] << " " << solide[it].discrete_gradient.col1[1] << " " << solide[it].discrete_gradient.col1[2] << endl;
	vtk << solide[it].discrete_gradient.col2[0] << " " << solide[it].discrete_gradient.col2[1] << " " << solide[it].discrete_gradient.col2[2] << endl;
	vtk << solide[it].discrete_gradient.col3[0] << " " << solide[it].discrete_gradient.col3[1] << " " << solide[it].discrete_gradient.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Epsilon_p
    vtk << "TENSORS epsilon_p double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++)
      {
	vtk << solide[it].epsilon_p.col1[0] << " " << solide[it].epsilon_p.col1[1] << " " << solide[it].epsilon_p.col1[2] << endl;
	vtk << solide[it].epsilon_p.col2[0] << " " << solide[it].epsilon_p.col2[1] << " " << solide[it].epsilon_p.col2[2] << endl;
	vtk << solide[it].epsilon_p.col3[0] << " " << solide[it].epsilon_p.col3[1] << " " << solide[it].epsilon_p.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Deformation plastique cumulée
    vtk << "SCALARS p double 1" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++)
      {
	vtk << solide[it].def_plas_cumulee << endl;
      }
    }
    vtk << "\n";
    //Rotation en x
    /*vtk << "VECTORS e double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++)
      {
	vtk << solide[it].e[0] << " " << solide[it].e[1] << " " << solide[it].e[2] << endl;
      }
    }
    vtk << "\n";*/
    //Vitesse de rotation
    /*vtk << "VECTORS omega double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(int it=0; it<nb_part; it++){
      for(int l= 0; l<solide[it].triangles.size(); l++)
      {
	vtk << solide[it].omega[0] << " " << solide[it].omega[1] << " " << solide[it].omega[2] << endl;
      }
    }
    vtk << "\n";*/
    vtk.close();
  }
  //Sortie sans reconstruction
  else{
    std::ostringstream oss;
    oss << "resultats/solide" << n << ".vtk";
    string s = oss.str();
    //cout << s << endl;
    const char* const solidevtk = s.c_str();
    //}
    
    //Ouverture des flux en donne en ecriture
    std::ofstream vtk;
    vtk.open(solidevtk,ios::out);
    if(vtk.is_open())
    {
      // cout <<"ouverture de xt.vtk reussie" << endl;
    } else {
      cout <<"ouverture de solide" << n << ".vtk rate" << endl;
    }
    //vtk << setprecision(15);
    int nb_points = 0;
    int nb_faces = 0;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      nb_faces += (*P).faces.size();
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	nb_points += (*F).vertex.size();
      }
    }
    
    //Initialisation du fichier vtk
    vtk << "# vtk DataFile Version 3.0" << endl;
    vtk << "#Simulation Euler" << endl;
    vtk << "ASCII" << endl;
    vtk<<"\n";
    vtk << "DATASET UNSTRUCTURED_GRID" << endl;
    vtk << "POINTS " << nb_points << " DOUBLE" << endl;
    
    //Sortie des points
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	for(std::vector<Vertex>::iterator V=(*F).vertex.begin();V!=(*F).vertex.end();V++){
	  vtk << (*P).mvt_t((*V).pos) << endl;
	}
      }
    }
    vtk << "\n";
    //Sortie des faces
    int point_tmp=0;
    vtk << "CELLS " << nb_faces << " " << nb_points+nb_faces << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*F).vertex.size();
	for(std::vector<Vertex>::iterator V=(*F).vertex.begin();V!=(*F).vertex.end();V++){
	  vtk << " " << point_tmp;
	  point_tmp++;
	}
	vtk << endl;
      }
    }
    vtk << "\n";
    vtk << "CELL_TYPES " << nb_faces << endl;
    for(int i=0;i<nb_faces;i++){
      vtk << 7 << endl;
    }
    vtk << "\n";
    vtk << "CELL_DATA " << nb_faces << endl;
    //Deplacement
    vtk << "VECTORS deplacement double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).Dx << endl;
      }
    }
    vtk << "\n";
    //Vitesse
    vtk << "VECTORS vitesse double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).u << endl;
      }
    }
    vtk << "\n";
    //Contrainte
    vtk << "TENSORS contraintes double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).contrainte.col1[0] << " " << (*P).contrainte.col1[1] << " " << (*P).contrainte.col1[2] << endl;
	vtk << (*P).contrainte.col2[0] << " " << (*P).contrainte.col2[1] << " " << (*P).contrainte.col2[2] << endl;
	vtk << (*P).contrainte.col3[0] << " " << (*P).contrainte.col3[1] << " " << (*P).contrainte.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Déformations
    vtk << "TENSORS deformations double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).discrete_gradient.col1[0] << " " << (*P).discrete_gradient.col1[1] << " " << (*P).discrete_gradient.col1[2] << endl;
	vtk << (*P).discrete_gradient.col2[0] << " " << (*P).discrete_gradient.col2[1] << " " << (*P).discrete_gradient.col2[2] << endl;
	vtk << (*P).discrete_gradient.col3[0] << " " << (*P).discrete_gradient.col3[1] << " " << (*P).discrete_gradient.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Epsilon_p
    vtk << "TENSORS epsilon_p double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).epsilon_p.col1[0] << " " << (*P).epsilon_p.col1[1] << " " << (*P).epsilon_p.col1[2] << endl;
	vtk << (*P).epsilon_p.col2[0] << " " << (*P).epsilon_p.col2[1] << " " << (*P).epsilon_p.col2[2] << endl;
	vtk << (*P).epsilon_p.col3[0] << " " << (*P).epsilon_p.col3[1] << " " << (*P).epsilon_p.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Deformation plastique cumulée
    vtk << "SCALARS p double 1" << endl;
    vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).def_plas_cumulee << endl;
      }
    }
    vtk << "\n";
    //Rotation en x
    /*vtk << "VECTORS e double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).e << endl;
      }
    }
    vtk << "\n";*/
    //Rotation en x
    /*vtk << "VECTORS eref double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).eref << endl;
      }
    }
    vtk << "\n";*/
    //Vitesse de rotation
    /*vtk << "VECTORS omega double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::vector<Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(*P).faces.begin();F!=(*P).faces.end();F++){
	vtk << (*P).omega << endl;
      }
    }
    vtk << "\n";*/
    vtk.close();
  }
}

// void Solide::breaking_criterion(){
// 	for(int it=0; it<solide.size(); it++){
// 		//cout<<"particule it "<<it<<endl; 
// 		for(int i=0; i<solide[it].faces.size(); i++){
// 			if(solide[it].faces[i].voisin >= 0){
// 				for(int iter=0; iter<solide.size(); iter++){
// 						if(it!=iter){
// 								if(solide[it].faces[i].voisin == iter){
// 									double distance = sqrt((CGAL::squared_distance(Point_3(solide[it].Dx.x() + solide[it].x0.x(), solide[it].Dx.y() + solide[it].x0.y(), solide[it].Dx.z() + solide[it].x0.z())  ,Point_3(solide[iter].Dx.x() + solide[iter].x0.x(),solide[iter].Dx.y() + solide[iter].x0.y(), solide[iter].Dx.z() + solide[iter].x0.z()) )));
// 	                  if( (distance - solide[it].faces[i].D0)/solide[it].faces[i].D0 >= k_max){
// 										cout<<"BREAK!!!!"<<endl; //cout<<"particule iter "<<iter<<endl;
// 										solide[it].faces[i].voisin = -2;
// 										int j;
// 										for(int f=0; f<solide[iter].faces.size(); f++){ //
// 											if(solide[iter].faces[f].voisin == it){
// 										    solide[iter].faces[f].voisin = -2;
// 												j=f;
// 											}
// 										}
// 										for(int count=0; count<solide[it].faces.size() ; count++){
// 											for(int ii=0; ii<solide[it].faces[count].vertex.size(); ii++)
// 											{ 
// 												std::vector<int> particules;
// 												for(int part=0; part<solide[it].faces[count].vertex[ii].particules.size(); part++){
// 													if(solide[it].faces[count].vertex[ii].particules[part] != iter){
// 														particules.push_back(solide[it].faces[count].vertex[ii].particules[part]);
// 													}
// 												}
// 													solide[it].faces[count].vertex[ii].particules.erase(solide[it].faces[count].vertex[ii].particules.begin(),
// 																																					solide[it].faces[count].vertex[ii].particules.end());
// 												  solide[it].faces[count].vertex[ii].particules = particules;
																																					
// 											}
// 									 }
// 									 for(int count=0; count<solide[iter].faces.size() ; count++){
// 										 for(int ii=0; ii<solide[iter].faces[count].vertex.size(); ii++)
// 											{ 
// 												std::vector<int> particules;
// 												for(int part=0; part<solide[iter].faces[count].vertex[ii].particules.size(); part++){
// 													if(solide[iter].faces[count].vertex[ii].particules[part] != it){
// 														particules.push_back(solide[iter].faces[count].vertex[ii].particules[part]);
// 													}
// 												}
// 												solide[iter].faces[count].vertex[ii].particules.erase(solide[iter].faces[count].vertex[ii].particules.begin(),
// 																																							solide[iter].faces[count].vertex[ii].particules.end());
// 												solide[iter].faces[count].vertex[ii].particules = particules;
// 											}
// 									  }
// 									  for(int count=0; count<solide.size() ; count++){
// 											if(count != it && count != iter){
// 												bool voisin_it = false;
// 												for(int jj=0; jj<solide[count].faces.size() ; jj++){
// 													if(solide[count].faces[jj].voisin==it){
// 														voisin_it=true;
// 													}
// 												}
// 												if(!voisin_it){
// 													//cout<<"voisin it"<<endl;
// 													for(int kk=0; kk<solide[it].faces[i].size() ; kk++){
// 															for(int jj=0; jj<solide[count].faces.size() ; jj++){
// 																for(int ii=0; ii<solide[count].faces[jj].vertex.size(); ii++)
// 																{   
// 																	if(solide[it].faces[i].vertex[kk].num == solide[count].faces[jj].vertex[ii].num){
// 																			std::vector<int> particules;
// 																			for(int part=0; part<solide[count].faces[jj].vertex[ii].particules.size(); part++){
// 																				if(solide[count].faces[jj].vertex[ii].particules[part] != it){
// 																					particules.push_back(solide[count].faces[jj].vertex[ii].particules[part]);
// 																				}
// 																			}
// 																			solide[count].faces[jj].vertex[ii].particules.erase(solide[count].faces[jj].vertex[ii].particules.begin(), solide[count].faces[jj].vertex[ii].particules.end());
// 																			solide[count].faces[jj].vertex[ii].particules = particules;
// 																  }
// 																}
// 															}
// 													}
// 													for(int kk=0; kk<solide[it].faces[i].size() ; kk++){
// 														for(int jj=0; jj<solide[it].faces.size() ; jj++){
// 															for(int ii=0; ii<solide[it].faces[jj].vertex.size(); ii++)
// 															{   
// 																if(solide[it].faces[i].vertex[kk].num == solide[it].faces[jj].vertex[ii].num){
// 																	std::vector<int> particules;
// 																	for(int part=0; part<solide[it].faces[jj].vertex[ii].particules.size(); part++){
// 																		if(solide[it].faces[jj].vertex[ii].particules[part] != count){
// 																			particules.push_back(solide[it].faces[jj].vertex[ii].particules[part]);
// 																		}
// 																	}
// 																	solide[it].faces[jj].vertex[ii].particules.erase(solide[it].faces[jj].vertex[ii].particules.begin(), solide[it].faces[jj].vertex[ii].particules.end());
// 																	solide[it].faces[jj].vertex[ii].particules = particules;
// 																}
// 															}
// 														}
// 													}
// 												}
// 												bool voisin_iter = false;
// 												for(int jj=0; jj<solide[count].faces.size() ; jj++){
// 													if(solide[count].faces[jj].voisin == iter){
// 														voisin_iter=true;
// 													}
// 												}
// 												if(!voisin_iter){
// 													//cout<<"voisin iter"<<endl;
// 													for(int kk=0; kk<solide[iter].faces[j].size() ; kk++){
// 														for(int jj=0; jj<solide[count].faces.size() ; jj++){
// 															for(int ii=0; ii<solide[count].faces[jj].vertex.size(); ii++)
// 															{   
// 																if(solide[iter].faces[j].vertex[kk].num == solide[count].faces[jj].vertex[ii].num){
// 																	std::vector<int> particules;
// 																	for(int part=0; part<solide[count].faces[jj].vertex[ii].particules.size(); part++){
// 																		if(solide[count].faces[jj].vertex[ii].particules[part] != iter){
// 																			particules.push_back(solide[count].faces[jj].vertex[ii].particules[part]);
// 																		}
// 																	}
// 																	solide[count].faces[jj].vertex[ii].particules.erase(solide[count].faces[jj].vertex[ii].particules.begin(), solide[count].faces[jj].vertex[ii].particules.end());
// 																	solide[count].faces[jj].vertex[ii].particules = particules;
// 																}
// 															}
// 														}
// 													}
// 													for(int kk=0; kk<solide[iter].faces[j].size() ; kk++){
// 														for(int jj=0; jj<solide[iter].faces.size() ; jj++){
// 															for(int ii=0; ii<solide[iter].faces[jj].vertex.size(); ii++)
// 															{   
// 																if(solide[iter].faces[j].vertex[kk].num == solide[iter].faces[jj].vertex[ii].num){
// 																	std::vector<int> particules;
// 																	for(int part=0; part<solide[iter].faces[jj].vertex[ii].particules.size(); part++){
// 																		if(solide[iter].faces[jj].vertex[ii].particules[part] != count){
// 																			particules.push_back(solide[iter].faces[jj].vertex[ii].particules[part]);
// 																		}
// 																	}
// 																	solide[iter].faces[jj].vertex[ii].particules.erase(solide[iter].faces[jj].vertex[ii].particules.begin(), solide[iter].faces[jj].vertex[ii].particules.end());
// 																	solide[iter].faces[jj].vertex[ii].particules = particules;
// 																}
// 															}
// 														}
// 													}
// 												}
												
// 											} //end !it et !iter
// 										}//end count
									  
// 									} //break
// 								}
// 						}
// 					}
// 			}
// 		}
// 	}
// 	//cout<<"Affiche "<<endl; Affiche();
// 	//cout<<"fin break"<<endl;
// }


#endif
