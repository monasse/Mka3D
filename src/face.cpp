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

Face::Face() : I_Dx(), vec_tangent_1(), vec_tangent_2(), voisins()
{
  centre = Point_3(0.,0.,0.);
  normale = Vector_3(1.,0.,0.);
  S = 0.;
  type = 0;
  id = -1; //Face pas remplie
}

Face & Face:: operator=(const Face &F){
  assert(this != &F);
  centre = F.centre;
  normale = F.normale;
  type = F.type;
  vertex.resize(F.vertex.size());
  for(int i= 0; i<F.vertex.size(); i++){
    vertex[i] = F.vertex[i];
  }
  for(int i= 0; i<F.voisins.size(); i++){
    voisins[i] = F.voisins[i];
  }
  for(int i= 0; i<F.c_reconstruction.size(); i++){
    c_reconstruction[i] = F.c_reconstruction[i];
  }
}

bool operator==(const Face &F1, const Face &F2) { //Compare les faces
  if(F1.type == 2) { //Triangle
    if(F1.vertex[0] == F2.vertex[0] && F1.vertex[1] == F2.vertex[1]  && F1.vertex[2] == F2.vertex[2])
      return true;
    else if(F1.vertex[0] == F2.vertex[0] && F1.vertex[1] == F2.vertex[2]  && F1.vertex[2] == F2.vertex[1])
      return true;
    else if(F1.vertex[0] == F2.vertex[1] && F1.vertex[1] == F2.vertex[0]  && F1.vertex[2] == F2.vertex[2])
      return true;
    else if(F1.vertex[0] == F2.vertex[1] && F1.vertex[1] == F2.vertex[2]  && F1.vertex[2] == F2.vertex[0])
      return true;
    else if(F1.vertex[0] == F2.vertex[2] && F1.vertex[1] == F2.vertex[0]  && F1.vertex[2] == F2.vertex[1])
      return true;
    else if(F1.vertex[0] == F2.vertex[2] && F1.vertex[1] == F2.vertex[1]  && F1.vertex[2] == F2.vertex[0])
      return true;
    else
      return false;
  }
  else if(F1.type == 3) {//Quad. Reprendre cette version !
    if(F1.vertex[0] != F2.vertex[0] && F1.vertex[0] != F2.vertex[1] && F1.vertex[0] != F2.vertex[2] && F1.vertex[0] != F2.vertex[3])
      return false;
    else {
      if(F1.vertex[1] != F2.vertex[0] && F1.vertex[1] != F2.vertex[1] && F1.vertex[1] != F2.vertex[2] && F1.vertex[1] != F2.vertex[3])
	return false;
      else {
	if(F1.vertex[2] != F2.vertex[0] && F1.vertex[2] != F2.vertex[1] && F1.vertex[2] != F2.vertex[2] && F1.vertex[2] != F2.vertex[3])
	  return false;
	else
	  return true;
      }
    }
  }
}

void Face::comp_quantities(Solide* Sol) { //, const Point_3& ext) {
  const Point_3 v1 = Sol->vertex[vertex[0]].pos;
  const Point_3 v2 = Sol->vertex[vertex[1]].pos;
  const Point_3 v3 = Sol->vertex[vertex[2]].pos;
  std::vector<Point_3> aux;
  aux.push_back(v1);
  aux.push_back(v2);
  aux.push_back(v3);
  if(type == 3)
    aux.push_back(Sol->vertex[vertex[2]].pos);
  centre = centroid(aux.begin(),aux.end());
  if(type == 2)
    S = 1./2. * sqrt(cross_product(Vector_3(v1,v2),Vector_3(v1,v3)).squared_length());
  else if(type == 3)
    S = 1./2. * sqrt(cross_product(Vector_3(v1,v2),Vector_3(v1,v3)).squared_length()) + 1./2. * sqrt(cross_product(Vector_3(v1,v2),Vector_3(v1,Sol->vertex[vertex[2]].pos)).squared_length());
  normale = orthogonal_vector(aux[0], aux[1], aux[2]);
  double norm = sqrt((normale.squared_length()));
  normale = normale / norm;
  D0 = 1000000000.;

  /*if(BC != 0) { //Pour les faces au bord, on calcule les vecteurs tangents
    double eps = pow(10., -14.);//std::numeric_limits<double>::epsilon(); //1e-14
    //Choix initial d'un repere orthonorme de la face
    Vector_3 s,t;
    if(normale[0]!=0. || normale[1]!=0.){
      s = Vector_3(-normale[1],normale[0],0.);
      s = s/(sqrt((s.squared_length())));
      t = cross_product(normale,s);
    }
    else if(normale[1]!=0. || normale[2]!=0.){
      s = Vector_3(0.,-normale[2],normale[1]);
      s = s/(sqrt((s.squared_length())));
      t = cross_product(normale,s);
    }
    else if(normale[0]!=0. || normale[2]!=0.){
      s = Vector_3(-normale[2],0,normale[0]);
      s = s/(sqrt((s.squared_length())));
      t = cross_product(normale,s);
    }
    else {
      s = Vector_3(0.,0.,1.);
      s = s/(sqrt((s.squared_length())));
      t = cross_product(normale,s); //t nul pk ?
    }
    //cout << "s : " << s << " t : " << t << endl;
    //Calcul de la matrice d'inertie de la face dans les deux axes a l'origine centre
    double Tss = 0.;
    double Ttt = 0.;
    double Tst = 0.;
    for(int i=0; i<vertex.size() ;i++){
      int ip = (i+1)%(vertex.size());
      Vector_3 V1(centre,Sol->vertex[vertex[i]].pos); //Modifier ici pour calcul ?
      Vector_3 V2(centre,Sol->vertex[vertex[ip]].pos);
      double As = (V1*s);
      double At = (V1*t);
      double Bs = (V2*s);
      double Bt = (V2*t);
      //cout << "As : " << As << " At : " << At << endl;
      Tss += 1./12.*(As*As+As*Bs+Bs*Bs);
      Ttt += 1./12.*(At*At+At*Bt+Bt*Bt);
      Tst += 1./24.*(2.*As*At+As*Bt+At*Bs+2.*Bs*Bt);
    }
    //cout << "Tst : " << Tst << endl;
    //Calcul des moments d'inertie
    double Delta = pow(Tss-Ttt,2)+4.*Tst*Tst;
    double Is,It; //Ce sont les 2 interties liées aux 2 vecs tangents dans le repère d'inertie de la face
    Is = (Tss+Ttt+sqrt(Delta))/2.;
    It = (Tss+Ttt-sqrt(Delta))/2.;
    //cout << "Moments : " << Is << " " << It << endl;
    //Diagonalisation
    if(abs(Tss-Ttt)>eps){
      if(sqrt((Tss-Is)*(Tss-Is))>eps){
	Vector_3 stemp = -Tst*s+(Tss-Is)*t;
	vec_tangent_1 = stemp/(sqrt((stemp.squared_length())));
	vec_tangent_2 = cross_product(normale,s);
      } else {
	Vector_3 stemp = -Tst*t+(Ttt-Is)*s;
	vec_tangent_1 = stemp/(sqrt((stemp.squared_length())));
	vec_tangent_2 = cross_product(normale,s);
      }
    } else {
      //cout << "Tst : " << sqrt(Tst*Tst) << " eps : " << eps << endl;
      if(sqrt(Tst*Tst)>eps){
	//cout << "On passe là" << endl;
	Vector_3 stemp = s+t;
	Vector_3 ttemp = -s+t;
	vec_tangent_1 = stemp/(sqrt((stemp.squared_length())));
	vec_tangent_2 = stemp/(sqrt((ttemp.squared_length())));
      }
    }
    }*/
}

#endif
