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

//Il faudra reprendre le calcul de la matrice d'inertie avec le formalisme Voronoi !!!!!

/*void Face::Inertie(const Particule& part){
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
    Vector_3 v1(centre,(part.vertices)[vertex[i]]);
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
}*/

#endif
