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

#include <stdio.h>
#include <math.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include "vitesse.cpp"
#include "geometry.hpp"

using namespace std;



int main()
{
   //Recover the parameters of the simulation from param.dat
  std::ifstream param("param.dat");
  if(!param){
    cout << "opening of param.dat failed" << endl;
  }
  string s;
  int numrep1, N_dim1, nimp1, Nmax1,Nx1,Ny1,Nz1,BCXl1,BCXr1,BCYl1,BCYr1,BCZl1,BCZr1;
  double rho1,nu1,E1,T1,cfl1, Lx1,Ly1,Lz1, Amort, KIc;
  bool rep1, flag2d1, rec1;
  param >> s >> rep1 >> s >> numrep1 >> s >> N_dim1 >> s >> flag2d1 >> s >> rho1 >> s >> nu1 >> s >> E1 >> s >> KIc >> s >> T1 >> s >> cfl1 >> s >> nimp1 >> s >> Nmax1 >> s >> rec1 >> s >> Amort >> s >> Lx1 >> s >> Ly1 >> s >> Lz1 >> s >>Nx1 >> s >> Ny1 >> s >> Nz1 >> s >> BCXl1 >> BCXr1 >> s >> BCYl1 >> BCYr1 >>s >> BCZl1 >> BCZr1;
  const double Lx = Lx1;
  const double Ly = Ly1;
  const double Lz = Lz1;
  const int Nx = Nx1;
  const int Ny = Ny1;
  const int Nz = Nz1;
  const double dx = Lx/Nx;
  const double dy = Ly/Ny;
  const double dz = Lz/Nz;
  const int BCXl = BCXl1;
  const int BCXr = BCXr1;
  const int BCYl = BCYl1;
  const int BCYr = BCYr1;
  const int BCZl = BCZl1;
  const int BCZr = BCZr1;

  //Opening of the output mesh file
  ofstream maillage;
  maillage.open("maillage.dat",ios::out);
  if(!maillage.is_open()){
    throw std::invalid_argument( "maillage.dat is not open !" );
  }
  //Printing vertices
  maillage << "POINTS " << (Nx+1)*(Ny+1)*(Nz+1) << endl;
  for(int i=0;i<Nx+1;i++){
    for(int j=0;j<Ny+1;j++){
      for(int k=0;k<Nz+1;k++){
	maillage << i*dx << " " << j*dy << " " << k*dz << endl;
      }
    }
  }
  maillage << endl;
  //Printing particles
  maillage << "SOLIDE " << Nx*Ny*Nz << endl;
  for(int i=0;i<Nx;i++){
    for(int j=0;j<Ny;j++){
      for(int k=0;k<Nz;k++){
	maillage << endl;
	maillage << "PARTICULE " << 6 << " ";
	//Boundary conditions
	int BC = 0; //Ajouter des possibilités de BC !!!!
	if(i==0){
	  BC = BCXl;
	}
	if(i==Nx-1 && BC!=1){
	  BC = BCXr;
	}
	if(j==0 && BC!=1){
	  BC = BCYl;
	}
	if(j==Ny-1 && BC!=1){
	  BC= BCYr;
	}
	if(k==0 && BC!=1){
	  BC = BCZl;
	}
	if(k==Nz-1 && BC!=1){
	  BC = BCZr;
	}
	maillage << BC << endl;
	//Position of the center of mass
	Point_3 P(dx*(i+0.5),dy*(j+0.5),dz*(k+0.5));
	maillage << "POSITION " << P << endl;
	if(BC==0){
	  maillage << "VITESSE " << velocity(P) << endl;
	  maillage << "VITROT " << omega(P) << endl;
	} else {
	  maillage << "VITESSE " << Vector_3(0,0,0) << endl;
	  maillage << "VITROT " << Vector_3(0,0,0) << endl;
	}
	//Connectivity
	//Face below
	maillage << 4 << " " << (i*(Ny+1)+j)*(Nz+1)+k << " " << (i*(Ny+1)+j+1)*(Nz+1)+k << " " << ((i+1)*(Ny+1)+j+1)*(Nz+1)+k << " " << ((i+1)*(Ny+1)+j)*(Nz+1)+k << " ";
	if(k>0){
	  maillage << (i*Ny+j)*Nz+k-1 << endl;
	} else {
	  maillage << -1 << endl;
	}
	//Face in front
	maillage << 4 << " " << (i*(Ny+1)+j)*(Nz+1)+k << " " << ((i+1)*(Ny+1)+j)*(Nz+1)+k << " " << ((i+1)*(Ny+1)+j)*(Nz+1)+k+1 << " " << (i*(Ny+1)+j)*(Nz+1)+k+1 << " ";
	if(j>0){
	  maillage << (i*Ny+j-1)*Nz+k << endl;
	} else {
	  maillage << -1 << endl;
	}
	//Face on the left
	maillage << 4 << " " << (i*(Ny+1)+j)*(Nz+1)+k << " " << (i*(Ny+1)+j)*(Nz+1)+k+1 << " " << (i*(Ny+1)+j+1)*(Nz+1)+k+1 << " " << (i*(Ny+1)+j+1)*(Nz+1)+k << " ";
	if(i>0){
	  maillage << ((i-1)*Ny+j)*Nz+k << endl;
	} else {
	  maillage << -1 << endl;
	}
	//Face on the right
	maillage << 4 << " " << ((i+1)*(Ny+1)+j)*(Nz+1)+k << " " << ((i+1)*(Ny+1)+j+1)*(Nz+1)+k << " " << ((i+1)*(Ny+1)+j+1)*(Nz+1)+k+1 << " " << ((i+1)*(Ny+1)+j)*(Nz+1)+k+1 << " ";
	if(i<Nx-1){
	  maillage << ((i+1)*Ny+j)*Nz+k << endl;
	} else {
	  maillage << -1 << endl;
	}
	//Face behind
	maillage << 4 << " " << (i*(Ny+1)+j+1)*(Nz+1)+k << " " << (i*(Ny+1)+j+1)*(Nz+1)+k+1 << " " << ((i+1)*(Ny+1)+j+1)*(Nz+1)+k+1 << " " << ((i+1)*(Ny+1)+j+1)*(Nz+1)+k << " ";
	if(j<Ny-1){
	  maillage << (i*Ny+j+1)*Nz+k << endl;
	} else {
	  maillage << -1 << endl;
	}
	//Face above
	maillage << 4 << " " << (i*(Ny+1)+j)*(Nz+1)+k+1 << " " << ((i+1)*(Ny+1)+j)*(Nz+1)+k+1 << " " << ((i+1)*(Ny+1)+j+1)*(Nz+1)+k+1 << " " << (i*(Ny+1)+j+1)*(Nz+1)+k+1 << " ";
	if(k<Nz-1){
	  maillage << (i*Ny+j)*Nz+k+1 << endl;
	} else {
	  maillage << -1 << endl;
	}
      }
    }
  }
  
	
  return 0;
  
 	
}
