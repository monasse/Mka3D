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
#include "particule.hpp"
#include "face.hpp"
#include "vertex.hpp"
#include "geometry.hpp"
#include "vitesse.hpp"
#include "forces_ext.hpp"
#include <iostream>
#include <map>
#include <set>
#include <string>
#ifndef SOLIDE_CPP
#define SOLIDE_CPP

Solide::Solide(const double& E, const double& nu){
  lambda = E * nu / (1.+nu) / (1. - 2.*nu);
  mu = E / 2. / (1.+nu);
}

Solide::Solide(){
  lambda = 0.;
  mu = 0.;
}

Solide::~Solide(){   
}

Solide & Solide::operator=(const Solide &S){
  assert(this != &S);
  //solide.resize(S.solide.size());
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    solide[P->first]= P->second;
  }
}

void Solide::Init(const char* s1, const bool& rep, const int& numrep, const double& rho){
  std::ifstream maillage(s1,ios::in);
  if(not(maillage))
    cout <<"ouverture du maillage ratee" << endl;

  //Importation des vertex
  string s;
  char aux;
  string ligne;
  while(getline(maillage, ligne) && ligne != "$Nodes") {} //On fait rien...
  getline(maillage, ligne); //Nombre de vertex. Pas utile.
  while(getline(maillage, ligne) && ligne != "$EndNodes") {
    istringstream  stm(ligne);
    int id;
    double x,y,z;
    stm >> id >> x >> y >> z;
    vertex.push_back(Vertex(Point_3(x,y,z))); //Vertex sont donnés dans l'ordre
  }

  //Importation des Particules
  getline(maillage, ligne); //$Elements
  getline(maillage, ligne); //Nbr Elements
  while(getline(maillage, ligne) && ligne != "$EndElements") {
    istringstream  stm(ligne);
    int id,type,nbr_tag,tag_1,tag_2;
    stm >> id >> type;
    if(type == 4) {
      int v1,v2,v3,v4;
      stm >> nbr_tag >> tag_1 >> tag_2 >> v1 >> v2 >> v3 >> v4;
      //Ajout des vertex de la particule
      Particule p(id);
      p.vertices.push_back(v1);
      p.vertices.push_back(v2);
      p.vertices.push_back(v3);
      p.vertices.push_back(v4);

      //Calcul des quantités volumiques (liées particule)
      p.barycentre(); //Calcul du barycentre
      p.volume(); //calcul du volume
      p.m = rho * p.V;
      //Comment fixer BC ?

      //Création des faces
      //face 1
      Face face1();
      face1.vertex.push_back(v2);
      face1.vertex.push_back(v3);
      face1.vertex.push_back(v4);
      std::set<Face>::iterator it = faces.find(face1);
      if(it == faces.end()) { //Face pas encore dans le set
	face1.comp_normal(v1); //Calcul de la normale sortante
	face1.surf(vertex[v2], vertex[v3], vertex[v4]); //Calcul surface face et centre
	face1.parts.push_back(id); //Ajout de la particule dans la face
	p.faces.push_back(&face1);
      }
      else {
        (it->parts).push_back(id); //Ajout de cette particule aux côtés de l'autre
	p.faces.push_back(it);
	Vector_3 aux(solide[(it->parts)[0]], x0);
	it->D0 = sqrt(aux.squared_length());
      }

      //face 2
      Face face2();
      face2.vertex.push_back(v1);
      face2.vertex.push_back(v3);
      face2.vertex.push_back(v4);
      std::set<Face>::iterator it = faces.find(face2);
      if(it == faces.end()) { //Face pas encore dans le set
	face2.comp_normal(v2); //Calcul de la normale sortante
	face2.surf(); //Calcul surface face
	face2.parts.push_back(id); //Ajout de la particule dans la face
	p.faces.push_back(&face2); //Ajout de la face dans la particule
      }
      else {
        (it->parts).push_back(id); //Ajout de cette particule aux côtés de l'autre
	p.faces.push_back(it);
	Vector_3 aux(solide[(it->parts)[0]], x0);
	it->D0 = sqrt(aux.squared_length());
      }

      //face 3
      Face face3();
      face3.vertex.push_back(v1);
      face3.vertex.push_back(v2);
      face3.vertex.push_back(v4);
      std::set<Face>::iterator it = faces.find(face3);
      if(it == faces.end()) { //Face pas encore dans le set
	face3.comp_normal(v3); //Calcul de la normale sortante
	face3.surf(); //Calcul surface face
	face3.parts.push_back(id); //Ajout de la particule dans la face
	p.faces.push_back(&face3);
      }
      else {
        (it->parts).push_back(id); //Ajout de cette particule aux côtés de l'autre
	p.faces.push_back(it);
	Vector_3 aux(solide[(it->parts)[0]], x0);
	it->D0 = sqrt(aux.squared_length());
      }

      //face 4
      Face face4();
      face4.vertex.push_back(v1);
      face4.vertex.push_back(v2);
      face4.vertex.push_back(v3);
      std::set<Face>::iterator it = faces.find(face4);
      if(it == faces.end()) { //Face pas encore dans le set
	face4.comp_normal(v4); //Calcul de la normale sortante
	face4.surf(); //Calcul surface face
	face4.parts.push_back(id); //Ajout de la particule dans la face
	p.faces.push_back(&face4);
      }
      else {
        (it->parts).push_back(id); //Ajout de cette particule aux côtés de l'autre
	p.faces.push_back(it);
	Vector_3 aux(solide[(it->parts)[0]], x0);
	it->D0 = sqrt(aux.squared_length());
      }   

      //Ajout de la particule dans le solide
      solide[id] = p;
    }
  }
}

void Solide::Solve_position(const double& dt, const bool& flag_2d, const double& t, const double& T){
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).solve_position(dt, flag_2d, t, T);
  }
  //breaking_criterion();
  /*update_triangles();
  for(int i=0;i<size();i++){
    for(std::vector<Triangle_3>::iterator it=solide[i].triangles.begin();it!=solide[i].triangles.end();it++){
      for(int k=0;k<3;k++){
	solide[i].bbox = Bbox(min(solide[i].bbox.xmin(),((*it).vertex(k).x())),min(solide[i].bbox.ymin(),((*it).vertex(k).y())),min(solide[i].bbox.zmin(),((*it).vertex(k).z())),max(solide[i].bbox.xmax(),((*it).vertex(k).x())),max(solide[i].bbox.ymax(),((*it).vertex(k).y())),max(solide[i].bbox.zmax(),((*it).vertex(k).z())));
      }
    }
    }*/	
}

void Solide::Solve_vitesse(const double& dt, const bool& flag_2d, const double& Amort, const double& t, const double& T){
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).solve_vitesse(dt, flag_2d, Amort, t , T);
  }
}

void Solide::Forces(const int& N_dim, const double& nu, const double& E, const double& dt, const double& t, const double& T){
  Forces_internes(N_dim,nu,E, dt);
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
      (P->second).Fi = (P->second).Fi + Forces_externes(t,T, *F);
  }
}

void Solide::Forces_internes(const int& N_dim, const double& nu, const double& E, const double& dt){
  bool plastifie = false;
  
  //Calcul de la contrainte dans chaque particule
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    //(P->second).Volume_libre();
    (P->second).discrete_gradient.col1 = Vector_3(0., 0., 0.); //Remet tous les coeffs de la matrice à 0.
    (P->second).discrete_gradient.col2 = Vector_3(0., 0., 0.);
    (P->second).discrete_gradient.col3 = Vector_3(0., 0., 0.);
    for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++){
      if(F->voisin>=0){
	int part = F->voisin;
	Vector_3 nIJ = (*F).normale;
	Matrix Dij_n(tens_sym(solide[part].Dx + solide[part].u * dt/2. - (P->second).Dx - (P->second).u * dt/2.,  nIJ) ); //Quadrature au point milieu pour calcul des forces !
	(P->second).discrete_gradient += (*F).S / 2. * Dij_n / (P->second).V;

      }
    }
    //cout << "Trace dev Def : " << (((P->second).discrete_gradient).dev()).tr() << endl;
    (P->second).contrainte = lambda * ((P->second).discrete_gradient - (P->second).epsilon_p).tr() * unit() + 2*mu * ((P->second).discrete_gradient - (P->second).epsilon_p);
    //cout << "Trace dev Contrainte : " << (((P->second).contrainte).dev()).tr() << endl;

    //Mettre ces valeurs dans le param.dat !!!!!
    double B = 292000000.; //En Pa. JC.
    double n = .31; //JC.
    double A = 90000000.; //En Pa. Vient de JC
    double H = 0.;
	
    (P->second).seuil_elas = A; // + B * pow((P->second).def_plas_cumulee, n);

    if(((P->second).contrainte - H * (P->second).epsilon_p).VM() > (P->second).seuil_elas) { //On sort du domaine élastique.
      plastifie = true;
      //Matrix n_elas(((*P).contrainte).dev() / (((*P).contrainte).dev()).norme() ); //Normale au domaine élastique de Von Mises
      Matrix n_elas( 1. / (((P->second).contrainte).dev()).norme() * ((P->second).contrainte).dev() ); //Normale au domaine élastique de Von Mises
      /*if((*P).n_elas_prev == -n_elas)
	cout << "Chargement dans sens oppose !" << endl;
      (*P).n_elas_prev = n_elas;*/
      //cout << "Trace n_elas : " << n_elas.tr() << endl;
      //cout << "Norme n_elas : " << n_elas.norme() << endl;
      //double delta_p = pow(((*P).contrainte.VM() - A) / B, 1./n) - (*P).def_plas_cumulee;
      double delta_p = (((P->second).contrainte - H * (P->second).epsilon_p).VM() - A) / (2*mu + H);
      //(*P).def_plas_cumulee = pow(((*P).contrainte.VM() - A) / B, 1./n); //Nouvelle déformation plastique.
      //cout << "Def plastique cumulee : " << (*P).def_plas_cumulee << endl;
      (P->second).epsilon_p += delta_p * n_elas;
      //cout << "Trace def plas : " << ((P->second).epsilon_p).tr() << endl; //Pb ! Non-nulle !!!!
      //cout << "Norme def plas : " << ((P->second).epsilon_p).norme() << endl;
      
      //((P->second).contrainte.VM() - A) / (2*mu) * n_elas;  //* signe( (P->second).contrainte ); //Plasticité parfaite
      //(P->second).contrainte = A * signe( (P->second).contrainte );
    }
  }

  if(plastifie)
    cout << "Plastification dans ce pas de temps !" << endl;

  
  //Calcul des forces pour chaque particule
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).Fi = Vector_3(0.,0.,0.);
    for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++){
      if((*F).voisin>=0){
	int part = (*F).voisin;
	Vector_3 nIJ = (*F).normale;
        Vector_3 Fij_elas( (*F).S / 2. * ( ((P->second).contrainte + solide[part].contrainte).tr() / 2. * nIJ + ((P->second).contrainte + solide[part].contrainte) / 2. * nIJ ) ); //Force du lien IJ !
	//cout << "Force : " << Fij_elas << endl;

	(P->second).Fi = (P->second).Fi + Fij_elas; // * nIJ; //Force sur particule
	
      }
      //cout << "Force interne : " << (P->second).Fi << endl;
    }
  }
}

double Solide::Energie(const int& N_dim, const double& nu, const double& E){
  return Energie_cinetique()+Energie_potentielle(N_dim, nu, E);
}

double Solide::Energie_cinetique(){
  double E = 0.;
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    double u2 = ((P->second).u.squared_length());
    E += 1./2. * (P->second).m * u2;
  }
  //cout << "Energie cinetique : " << E << endl;
  return E;
}

double Solide::Energie_potentielle(const int& N_dim, const double& nu, const double& E){
  double Ep = 0.;

  double B = 292000000.; //En Pa. JC.
  double n = .31; //JC.
  double A = 90000000.; //En Pa. Vient de JC

  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    Ep += ( 0.5 * contraction_double(((P->second)).contrainte, ((P->second)).discrete_gradient - ((P->second)).epsilon_p)  + B * pow(((P->second)).def_plas_cumulee, 1. + n) / (n + 1.) + A * ((P->second)).def_plas_cumulee ) * ((P->second)).V;
  }
  //cout << "Energie potentielle : " << Ep << endl;
  return Ep;
}

double Solide::pas_temps(const double& t, const double& T, const double& cfls, const double& E, const double& nu, const double& rhos){
  double eps = 1e-14;//std::numeric_limits<double>::epsilon();
  double dt = std::numeric_limits<double>::infinity();
  //Restriction CFL liee aux forces internes
  double cs = sqrt(E*(1.-nu)/rhos/(1.+nu)/(1.-2.*nu));
  //Calcul du rayon de la sphï¿½re inscrite
  double sigma = 100000.;
  for(int i=0;i<size();i++){
    for(int j=0;j<solide[i].faces.size();j++){
      sigma = min(sigma,solide[i].faces[j].D0);
    }
  }
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    for(int j=0;j<(P->second).faces.size();j++){
      if((P->second).faces[j].voisin>=0){
	dt = min(dt,cfls*(P->second).faces[j].D0/cs);
	/*double Imin = min(min((P->second).I[0],(P->second).I[1]),(P->second).I[2]);
	double S = (P->second).faces[j].S;
	double D0 = (P->second).faces[j].D0;
	double kappa = 1.;
	double alphan = (2.+2.*nu-kappa)*E/4./(1.+nu)/S*((P->second).faces[j].Is+(P->second).faces[j].It);
	double alphas = E/4./(1.+nu)/S*((2.+2.*nu+kappa)*(P->second).faces[j].Is-(2.+2.*nu-kappa)*(P->second).faces[j].It);
	double alphat = E/4./(1.+nu)/S*((2.+2.*nu+kappa)*(P->second).faces[j].It-(2.+2.*nu-kappa)*(P->second).faces[j].Is);
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphan));
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphas));
	dt = min(dt,cfls*sqrt(Imin*D0/S/alphat));*/
      }
    }
  }
  dt = min(dt,T-t);
  return dt;
}

/*void Solide::update_triangles(){
  for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
    (P->second).triangles_prev = (P->second).triangles;
    (P->second).normales_prev = (P->second).normales;
    (P->second).fluide_prev = (P->second).fluide;
    for(int it=0;it<(P->second).triangles.size();it++){
      (P->second).Points_interface_prev[it] = (P->second).Points_interface[it];
      (P->second).Triangles_interface_prev[it] = (P->second).Triangles_interface[it];
      (P->second).Position_Triangles_interface_prev[it] = (P->second).Position_Triangles_interface[it];
      (P->second).Points_interface[it].erase((P->second).Points_interface[it].begin(),(P->second).Points_interface[it].end());
      (P->second).Triangles_interface[it].erase((P->second).Triangles_interface[it].begin(),(P->second).Triangles_interface[it].end());	(P->second).Position_Triangles_interface[it].erase((P->second).Position_Triangles_interface[it].begin(), (P->second).Position_Triangles_interface[it].end());
    }
    (P->second).triangles.erase((P->second).triangles.begin(),(P->second).triangles.end());
    (P->second).normales.erase((P->second).normales.begin(),(P->second).normales.end());
    (P->second).fluide.erase((P->second).fluide.begin(),(P->second).fluide.end());
    (P->second).vide.erase((P->second).vide.begin(),(P->second).vide.end());
		
    //Calcul de la nouvelle position des triangles
    for(int f=0;f<(P->second).faces.size();f++){
      Point_3 s,r,v,t;
			
      if((P->second).faces[f].size() == 3){
	vector<Point_3> ri,vi,si ;
	for(int part=0; part<(P->second).faces[f].vertex[0].size();part++){
	  int p = (P->second).faces[f].vertex[0].particules[part];
	  ri.push_back(solide[p].mvt_t((P->second).faces[f].vertex[0].pos));
	}
	r = centroid(ri.begin(),ri.end());
				
				
	for(int part=0;part<(P->second).faces[f].vertex[1].size();part++){
	  int p = (P->second).faces[f].vertex[1].particules[part];
	  vi.push_back(solide[p].mvt_t((P->second).faces[f].vertex[1].pos));
	}
	v = centroid(vi.begin(),vi.end());

	for(int part=0;part<(P->second).faces[f].vertex[2].size();part++){
	  int p = (P->second).faces[f].vertex[2].particules[part];
	  si.push_back(solide[p].mvt_t((P->second).faces[f].vertex[2].pos));
	}
	s = centroid(si.begin(),si.end());
				
	Vector_3 vect0(r,v);
	Vector_3 vect1(r,s);
	Triangle_3 Tri(r,v,s);
	(P->second).triangles.push_back(Tri);
	Vector_3 normale = cross_product(vect0,vect1);
	normale = normale*(1./sqrt((normale.squared_length())));
	(P->second).normales.push_back(normale);
	if((P->second).faces[f].voisin < 0){
	  (P->second).fluide.push_back(true);
	} else {
	  (P->second).fluide.push_back(false);
	}
	if( (P->second).faces[f].voisin == -2){
	  (P->second).vide.push_back(true);
	} else {
	  (P->second).vide.push_back(false);
	}
      }
      else{
			
	vector<Point_3> si;
	si.push_back((P->second).mvt_t((P->second).faces[f].centre));
	int j = (P->second).faces[f].voisin;
	if(j>=0){
	  si.push_back(solide[j].mvt_t((P->second).faces[f].centre));
	}
	s = centroid(si.begin(),si.end());
			
	for(int k=0;k<(P->second).faces[f].size();k++){
	  int kp = (k+1)%((P->second).faces[f].size());
	  vector<Point_3> ri,vi;
	  for(int part=0;part<(P->second).faces[f].vertex[k].size();part++){
	    int p = (P->second).faces[f].vertex[k].particules[part];
	    ri.push_back(solide[p].mvt_t((P->second).faces[f].vertex[k].pos));
	  }
	  r = centroid(ri.begin(),ri.end());
	  for(int part=0;part<(P->second).faces[f].vertex[kp].size();part++){
	    int p = (P->second).faces[f].vertex[kp].particules[part];
	    vi.push_back(solide[p].mvt_t((P->second).faces[f].vertex[kp].pos));
	  }
	  v = centroid(vi.begin(),vi.end());
	  Vector_3 vect0(s,r);
	  Vector_3 vect1(s,v);
	  Triangle_3 Tri(s,r,v);
	  (P->second).triangles.push_back(Tri);
	  Vector_3 normale = cross_product(vect0,vect1);
	  normale = normale*(1./sqrt((normale.squared_length())));			  
	  (P->second).normales.push_back(normale);
	  if((P->second).faces[f].voisin < 0){
	    (P->second).fluide.push_back(true);
	  } 
	  else {
	    (P->second).fluide.push_back(false);
	  }
	  if( (P->second).faces[f].voisin == -2){
	    (P->second).vide.push_back(true);
	  } else {
	    (P->second).vide.push_back(false);
	  }
	}
      }
			
    }//Calcul de la nouvelle position des triangles
  }
  }*/

void Solide::Impression(const int &n, const bool &reconstruction){ //Sortie au format vtk
  int nb_part = solide.size();
//Version avec reconstruction
  /*if(reconstruction){
    int nb_triangles = 0.;
    for(int it=0; it<nb_part; it++){
      nb_triangles += solide[it].triangles.size();
    }

//const char* solidevtk;
//{
    std::ostringstream oss;
    oss << "solide" << n << ".vtk";
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
    vtk.close();
  }
  //Sortie sans reconstruction
  else{*/
    std::ostringstream oss;
    oss << "solide" << n << ".vtk";
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
    int size = 0;
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      nb_faces += (P->second).faces.size();
      nb_points += (P->second).vertices.size();
      for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
	size += F->nb_vertex; //Egale à nb_points au final ?
    }
    size += nb_faces;
    
    //Initialisation du fichier vtk
    vtk << "# vtk DataFile Version 3.0" << endl;
    vtk << "#Simulation Euler" << endl;
    vtk << "ASCII" << endl;
    vtk<< "\n";
    vtk << "DATASET UNSTRUCTURED_GRID" << endl;
    vtk << "POINTS " << nb_points << " DOUBLE" << endl;
    
    //Sortie des points
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Point_3>::iterator V=(P->second).vertices.begin();V!=(P->second).vertices.end();V++) {
	vtk << (P->second).mvt_t(*V) << endl;
      }
    }
    vtk << "\n";
    
    //Sortie des faces
    int point_tmp=0;
    vtk << "CELLS " << nb_faces << " " << size << endl;
    int compteur_vertex = 0;
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
	vtk << F->nb_vertex;
	for(int k=0 ; k<F->nb_vertex ; k++)
	  vtk << " " << compteur_vertex + (F->vertex)[k];
	vtk << endl;
      }
      compteur_vertex += (P->second).vertices.size();
      //vtk << endl;
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
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
	vtk << (P->second).Dx << endl;
    }
    vtk << "\n";
    //Vitesse
    vtk << "VECTORS vitesse double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
	vtk << (P->second).u << endl;
    }
    vtk << "\n";
    //Contrainte
    vtk << "TENSORS contraintes double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      vtk << (P->second).contrainte.col1[0] << " " << (P->second).contrainte.col1[1] << " " << (P->second).contrainte.col1[2] << endl;
      vtk << (P->second).contrainte.col2[0] << " " << (P->second).contrainte.col2[1] << " " << (P->second).contrainte.col2[2] << endl;
      vtk << (P->second).contrainte.col3[0] << " " << (P->second).contrainte.col3[1] << " " << (P->second).contrainte.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Déformations
    vtk << "TENSORS deformations double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      vtk << (P->second).discrete_gradient.col1[0] << " " << (P->second).discrete_gradient.col1[1] << " " << (P->second).discrete_gradient.col1[2] << endl;
      vtk << (P->second).discrete_gradient.col2[0] << " " << (P->second).discrete_gradient.col2[1] << " " << (P->second).discrete_gradient.col2[2] << endl;
      vtk << (P->second).discrete_gradient.col3[0] << " " << (P->second).discrete_gradient.col3[1] << " " << (P->second).discrete_gradient.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Epsilon_p
    vtk << "TENSORS epsilon_p double" << endl;
    //vtk << "LOOKUP_TABLE default" << endl;
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++) {
      vtk << (P->second).epsilon_p.col1[0] << " " << (P->second).epsilon_p.col1[1] << " " << (P->second).epsilon_p.col1[2] << endl;
      vtk << (P->second).epsilon_p.col2[0] << " " << (P->second).epsilon_p.col2[1] << " " << (P->second).epsilon_p.col2[2] << endl;
      vtk << (P->second).epsilon_p.col3[0] << " " << (P->second).epsilon_p.col3[1] << " " << (P->second).epsilon_p.col3[2] << endl;
      }
    }
    vtk << "\n";
    //Deformation plastique cumulée
    vtk << "SCALARS p double 1" << endl;
    vtk << "LOOKUP_TABLE default" << endl;
    for(std::map<int, Particule>::iterator P=solide.begin();P!=solide.end();P++){
      for(std::vector<Face>::iterator F=(P->second).faces.begin();F!=(P->second).faces.end();F++)
	vtk << (P->second).def_plas_cumulee << endl;
    }
    vtk << "\n";
    vtk.close();
}

#endif
