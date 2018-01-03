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
#include <map>
#include <string>
#ifndef PARTICULE_CPP
#define PARTICULE_CPP

//using namespace std;

inline double signe(const double &x)
{
  return (x < 0.) ? -1. : 1. ;
}

//Function absolute value
double abs(const double &x)
{
  return max(x,-x);
}

Particule::Particule(const int& Id):discrete_gradient(), contrainte(), epsilon_p(), vertices() {
  id = Id;
  def_plas_cumulee = 0.; //Déformation plastique cumulée du lien
  seuil_elas = 0.;
  fixe = 0;
}

Particule::Particule():discrete_gradient(), contrainte(), epsilon_p()
{
  id = 0;
  def_plas_cumulee = 0.; //Déformation plastique cumulée du lien
  seuil_elas = 0.;
  fixe = 0;
}


Particule::~Particule(){
}

Particule & Particule:: operator=(const Particule &P){
  assert(this != &P);
  bbox = P.bbox;
  cube  = P.cube;
	
  faces = P.faces;
  fixe = P.fixe;
  m  = P.m; 
  V = P.V; 
  epsilon = P.epsilon; 
  /*for(int i=0; i<3;i++){
    I[i] = P.I[i];
    for(int j=0; j<3;j++){
    rotref[i][j] = P.rotref[i][j];
    }
    }*/
	
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
}
	
/*	triangles.resize(P.triangles.size()); 
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
	
	}*/

/*void Particule::Affiche(){
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

  }*/


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
  }


  //Dx = displacement_BC(x0, Dx, t, T);

  //Mise a jour de la transformation donnant le mouvement de la particule
  mvt_tprev = mvt_t;
  //Aff_transformation_3 rotation(rot[0][0],rot[1][0],rot[2][0],rot[0][1],rot[1][1],rot[2][1],rot[0][2],rot[1][2],rot[2][2]);
  Aff_transformation_3 translation(Vector_3(Point_3(0.,0.,0.),x0)+Dx);
  Aff_transformation_3 translation_inv(Vector_3(x0,Point_3(0.,0.,0.)));
  mvt_t = translation*(/*rotation*/translation_inv);
	//cout<<"position du centre de la particule "<<x0+Dx<<endl;
}

void Particule::solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  if(fixe==1){
    u = Vector_3(0.,0.,0.);
    //omega = Vector_3(0.,0.,0.);
  } else {
    if(fixe==0){
      u = u+(Fi+Ff)/2.*(dt/m)*Amort; // + velocity_BC(x0, t, T, Dx); //Conditions aux limites en vitesse ajoutées ici
    }
    else if(fixe==2 || fixe==3){
      u = velocity_BC(x0, t, T, Dx); //Vector_3(0.,0.,0.);
    }
  }

  /*if(x0.z() <= 4.)
    u = Vector_3(0,0,-10.); //En m.s^-1
  else if(x0.z() >= 14.)
    u= Vector_3(0,0,0);
  else
  u = u+(Fi+Ff)/2.*(dt/m)*Amort;*/
    
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
    Q[2][2] = 1.-2.*(eref[0]*eref[0]+eref[1]*eref[1]);
    double e0 = sqrt(abs(1.-(e.squared_length())));
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
  //}Fin du calcul dans le cas d'une particule libre
}

Vector_3 Particule::vitesse_parois(const Point_3& X_f){
		
  Vector_3 V_f = u_half + cross_product(omega_half, Vector_3(x0 + Dx,X_f));

  return V_f;
}	

Vector_3 Particule::vitesse_parois_prev(const Point_3& X_f){
	
  Vector_3 V_f = u_half + cross_product(omega_half, Vector_3(x0 + Dxprev,X_f));
	
  return V_f;
}	

/*void Face::compProjectionIntegrals(double &P1, double &Pa, double &Pb, double &Paa, double &Pab, double &Pbb, double &Paaa, double &Paab, double &Pabb, double &Pbbb, const int& a, const int& b, const int& c){
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
}*/

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



/*void Particule::Inertie(const double &rho){
  double eps = 1.e-14;//std::numeric_limits<double>::epsilon();
  double T1,Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Tzx;
  //CompVolumeIntegrals(T1,Tx,Ty,Tz,Txx,Tyy,Tzz,Txy,Tyz,Tzx);
  double R[3][3];
  double xG = (x0[0]);
  double yG = (x0[1]);
  double zG = (x0[2]);
  
  R[0][0] = rho*(Tyy-2.*yG*Ty+yG*yG*T1+Tzz-2.*zG*Tz+zG*zG*T1);
  R[1][0] = R[0][1] = rho*(Txy-yG*Tx-xG*Ty+xG*yG*T1);
  R[2][0] = R[0][2] = rho*(Tzx-zG*Tx-xG*Tz+xG*zG*T1);
  R[1][1] = rho*(Txx-2.*xG*Tx+xG*xG*T1+Tzz-2.*zG*Tz+zG*zG*T1);
  R[1][2] = R[2][1] = rho*(Tyz-zG*Ty-yG*Tz+yG*zG*T1);
  R[2][2] = rho*(Tyy-2.*yG*Ty+yG*yG*T1+Txx-2.*xG*Tx+xG*xG*T1);
  
  //Masse et volume
  //V = T1;
  //m = rho*T1;
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
  
  //Calcul des moments d'inertie des faces (pour le calcul des torsions)
  for(int i=0;i<faces.size();i++){
    faces[i].Inertie();
  }
}*/

#endif
