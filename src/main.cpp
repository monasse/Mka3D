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
 \mainpage
\authors Laurent Monasse and Adela Puscas 
 \section intro Introduction
 3d elastic solid simulation code (Mka3d)
The numerical method uses a Discrete Element method.

 \section description Description of the project
 
 To use the code, execute the following commands:
 - cmake .
 - make: compilation 
 - ./main: execution 

 A public repository for the code is on <a href="https://github.com/monasse/Mka3D">github</a>.

 Ce que vous devez remplir pour lancer une simulation:
 
 - \a parametres.hpp: d&eacute;finir les param&egrave;tres du probl&egrave;me; 
 - \a parametres.cpp: d&eacute;finir l'&eacute;tat initial du fluide: densit&eacute;, pression, vitesse;
 - cr&eacute;ation du fichier <b> "maillage.dat" </b>d&eacute;finissant le maillage pour le solide.\n
 
 
 
 
Cr&eacute;ation du fichier "maillage.dat": \n
- Premi&egrave;re ligne le mot cl&eacute; \b "POINTS" et le nombre de sommets du solide. \n
- Lignes suivantes : l'ensemble des coordonn&eacute;es de sommets. Chaque sommet sur une ligne.\n
- Le mot cl&eacute; \b "SOLIDE" et le nombre de particules solides du syst&egrave;me .\n
- Les caract&eacute;ristiques de chaque particule: \n
 - le mot cl&eacute; \b "PARTICULE" suivi du nombre de faces et la condition sur le mouvement de la particule: 0 ou 1. Si 0, la particule peut bouger librement si 1, la particule est fix&eacute;e &agrave; sa position initiale.\n
 - le mot cl&eacute; \b "POSITION" et les coordonn&eacute;es du centre de la particule. \n
 - le mot cl&eacute; \b "VITESSE" et les composantes de la vitesse initiale de la particule. \n
 - le mot cl&eacute; \b "VITROT"  et les composantes de la vitesse de rotation initiale de la particule. \n
- Description g&eacute;om&eacute;trique de chacune des particules face par face: 
 - nombre de sommets de la face, les num&eacute;ros de sommets (dans le m&ecirc;me ordre que dans le listing 
 des coordonn&eacute;es de sommets) d&eacute;crivant la face et un dernier nombre, part, indiquent le num&eacute;ro 
 de la particule (dans le m&ecirc;me ordre que dans le listing de particules) &eacute;ventuellement 
 pr&eacute;sente de l'autre c√¥t&eacute; de la face. Si aucune particule n'est pr&eacute;sente de l'autre c√¥t&eacute; 
 de la face (donc si la particule est en contact avec le fluide), on a part =-1. \n
 
\remark L'ordre de ces sommets donne l'orientation de la normale sortante &agrave; la particule,  
 les sommets sont num&eacute;rot&eacute;s de 0 &agrave; nombre de sommets-1 et les particules sont 
 num&eacute;rot&eacute;s de 0 &agrave; nombre de particules-1 comme dans un tableau C++. \n
 
 
Exemple: \n
 
			POINTS 8 \n
			0.4 0.4 0.4 \n
			0.9 0.4 0.4 \n
			0.9 0.6 0.4 \n
			0.4 0.6 0.4 \n
			0.4 0.4 0.6 \n
			0.9 0.4 0.6 \n
			0.9 0.6 0.6 \n
			0.4 0.6 0.6 \n
			
			SOLIDE 1 \n
			
			PARTICULE 6 0 \n
			POSITION  1 1 1 \n
			VITESSE 1. 0. 0.\n
			VITROT 0. 0. 0.\n
			4 0 3 2 1 -1 \n
			4 7 4 5 6 -1 \n
			4 0 4 7 3 -1 \n
			4 1 2 6 5 -1 \n
			4 4 0 1 5 -1 \n
			4 3 7 6 2 -1 \n
 
 
 Les r&eacute;sultats sont &eacute;crits dans le dossier \b resultats. Certains fichiers sont remplis &agrave; 
 chaque pas de temps :
 les fichiers energie.dat, solide_center.dat et temps_reprise.dat donnent respectivement l'&eacute;volution
 de l'&eacute;nergie, la position du centre du solide au cours du temps et le temps de la simulation. 
 Les fichiers fluide*.vtk et solide*.vtk sont &eacute;crits un nombre limit&eacute;
 de fois au cours du calcul; fluide*.vtk et solide*.vtk donnent l'&eacute;tat du fluide et la position du solide,
 peuvent &ecirc;tre lus avec Paraview. Le fichier temps.dat est &eacute;crit &agrave; la fin de la simulation 
 donnent le temps de calcul (temps machine) des differents sous-proc&eacute;dures  du code. \n
 Il est possible d'effectuer une reprise sur un calcul interrompu &agrave; partir des sorties fluide*.vtk
 et solide*.vtk: il suffit de changer le flag de reprise bool rep=false en bool rep=true dans le fichier parametres.h
 et d'indiquer &agrave; partir de quel point de reprise on souhaite reprendre &agrave; l'aide de la variable int numrep.
 Par exemple, pour reprendre avec fluide1.vtk et solide1.vtk (c'est-a-dire &agrave; partir de la premi&egrave;re 
 sortie), on utilisera int numrep=1.
 
 
 \remark Les proc&eacute;dures concernant le fluide se trouvent dans les fichiers fluide.hpp et fluide.cpp. Celles concernant le solide dans solide.hpp et solide.cpp. Les proc&eacute;dures r&eacute;alisant le couplage sont regroup&eacute;es dans les fichiers couplage.cpp, intersections.hpp et intersections.cpp. Les fichiers parametres.hpp et parametres.cpp sont d&eacute;dies aux d&eacute;finitions des param&eacute;tr&eacute;s du probl&egrave;me (param&eacute;tr&eacute;s physique, li&eacute;s aux maillages fluide et solide, au couplage, etc.). La r&eacute;solution du probl&egrave;me est r&eacute;alis&eacute;e  dans le fichier principal main.cpp. 
 
 
 */


/*!
*  \file main.cpp
*  \brief Fonction principale. Initialisation du probl&egrave;me et r&eacute;solution.
*/

#include <iostream>
#include <ctime>
#include "solide.hpp"
using namespace std;      

/*!
* \fn int main()
*\brief Initialisation du probl&egrave;me et r&eacute;solution:

- Initialisation du solide via la fonction \a Solide.Init(const char*) et du fluide via la fonction  Grille.Init().
- R&eacute;solution du probl&egrave;me:
 - R&eacute;solution des &eacute;quations fluides via la fonction \a Grille.Solve(const double, double, int).
 - Calcul des forces internes via la fonction \a Solide.Forces_internes().
 -  Calcul des forces (\a Particule.Ff) et moments fluides (\a Particule.Mf) exerc&eacute;s sur le solide via la fonction \a Grille.Forces_fluide(Solide&, const double).
 - Mise &agrave; jour de la position du solide via la fonction \a Solide.Solve_position(double).
 - Calcul de la vitesse du solide via la fonction \a Solide.Solve_vitesse(double).
 - Intersection de la grille fluide avec le solide via la fonction \a Grille.Parois(Solide&, double).
 - Calcul de la quantit&eacute; balay&eacute;e par le solide via la fonction \a Grille.Swap_2d(double, Solide&).
 - Modification des flux fluide via la fonction \a Grille.Modif_fnum(double).
 - M&eacute;lange conservatif de petites cellules coup&eacute;es via la fonction \a Grille.Mixage().
 - Remplissage des cellules fictives via la fonction \a Grille.Fill_cel(Solide&).
 - Imposition des conditions aux limites via la fonction \a Grille.BC().
 
*\return int
*/
int main(){

  //Recover the parameters of the simulation from param.dat
  std::ifstream param("param.dat");
  if(!param){
    cout << "opening of param.dat failed" << endl;
  }
  string s;
  int numrep1, N_dim1, nimp1, Nmax1;
  double rho1,nu1,E1,T1,cfl1,Amort;
  bool rep1, flag2d1, rec1;
  param >> s >> rep1 >> s >> numrep1 >> s >> N_dim1 >> s >> flag2d1 >> s >> rho1 >> s >> nu1 >> s >> E1 >> s >> T1 >> s >> cfl1 >> s >> nimp1 >> s >> Nmax1 >> s >> rec1 >> s >> Amort;
  const bool rep = rep1; //Recovery flag
  const int numrep = numrep1; //File number from which to possibly restart
  const int N_dim=N_dim1; //Number of dimensions of the problem
  const bool flag_2d = flag2d1; //Is the problem 2d ?
  const double rho = rho1;  //Density of the solid
  const double nu = nu1; //Poisson's ratio
  const double E = E1; //Young modulus
  const double T = T1;   //Total simulation time
  const double cfl = cfl1; //CFL condition number (must be less than 1)
  const int nimp = nimp1; //Number of outputs
  const double dtimp = T/nimp;        //Time-step between two consecutive outputs
  const int Nmax = Nmax1;           //Maximal number of time-steps
  const bool rec = rec1;   //Flag pour dÈcider si la sortie se fait avec reconstruction de l'interface ou particule par particule
  const double Amortissement = Amort; //Ajoute une force de frottement entre 0 et 1 pour amortir la solution
  
  char temps_it[]="temps.dat";
  char temps_reprise[]="temps_reprise.dat";

//En cas de reprise
  double temps[numrep+1];
  if(rep){
    //reprise();
    std::ifstream in(temps_reprise,ios::in);
    if(!in){
      cout <<"ouverture de temps_reprise.dat rate" << endl;
    }
    for(int i=0;i<numrep+1;i++){
      in >> temps[i];
    }
    //Recuperation de l'energie
    int result = system("cp energie.dat energie_reprise.dat");
    }
    std::ifstream in_energie("energie_reprise.dat",ios::in);
  
	
  //Ouverture des flux en donne en ecriture
  std::ofstream temps_iter(temps_it,ios::out);
  std::ofstream sorties_reprise(temps_reprise,ios::out);
  if(temps_iter)
  {
    // cout <<"ouverture de 'temps.dat' reussie" << endl;
  } else {
    cout <<"ouverture de 'temps.dat' rate" << endl;
  }
  if(rep){
    for(int i=0;i<numrep+1;i++){
      sorties_reprise << temps[i] << endl;
    }
  }
  
  
  char energie[]="energie.dat";
  
  //Ouverture des flux en donne en ecriture
  std::ofstream ener(energie,ios::out);
  if(ener)
  {
    // cout <<"ouverture de xt.vtk reussie" << endl;
  } else {
    cout <<"ouverture de energie.dat rate" << endl;
  }
  double dE0rep,dE0Srep,dm0;
  if(rep){
    double t_ener = 0.;
    double E,dE;
    for(int i=0;t_ener<temps[numrep];i++){
      in_energie >>t_ener >> E >>  dE;
      if(t_ener<=temps[numrep]){
	ener << t_ener << " " << E << " " << dE << endl;
	dE0rep = dE;
      }
    }
    int result = system("rm energie_reprise.dat");
  } 
  
  
  
  double t=0., dt=0.;
  if(rep){
    t = temps[numrep];
  }
  Solide S(E, nu);
  //Initialization from file "maillage*.dat", with possible restart depending on rep
  S.Init("packing.custom1", "packing.custom2", "packing.custom3", rep, numrep, rho);

  cout << "Lecture des fichiers de maillage terminÈe !" << endl;
  

  /*for(std::map<int, Particule>::iterator P=S.solide.begin();P!=S.solide.end();P++){
    cout << "Particle ID : " << (P->second).id << endl;
    }*/
  	
  //Initialization of time measurements
  int iter=0;	
  clock_t start,end;
  start =clock();
  
  int kimp = 0; //Output number
  double next_timp = dtimp; //Time of the next output
  if(rep){
    kimp = numrep;
    next_timp = t+dtimp;
  } else {
    S.Impression(kimp,rec);
    sorties_reprise << t << endl;
  }
  kimp++;
  
  double E0 = S.Energie(N_dim, nu, E);
  if(rep){
    E0 -= dE0rep;
  }
  S.Forces_internes(N_dim, nu, E, dt);
  int nb_part = S.size();

  //Iterations on the time-steps
  for (int n=0; (t<T) && n<Nmax; n++){
    cout<<"iteration="<<n<< " dt="<<dt<<" t="<<t<<endl;
    //Output at prescribed times
    if(t>next_timp){
      S.Impression(kimp,rec);
      sorties_reprise << t << endl;
      kimp++;
      next_timp += dtimp;
    }
    //Variation of energy
    cout<< "Solid energy:" << S.Energie(N_dim, nu, E) <<endl;
    //Variation of momentum
    Vector_3 qdm(0,0,0);
    for(std::map<int, Particule>::iterator P=S.solide.begin();P!=S.solide.end();P++){
      qdm = qdm+ (P->second).m*(P->second).u;
    }
    
    ener << t << " " << S.Energie(N_dim, nu, E) << " " << S.Energie(N_dim, nu, E)-E0 << " " << qdm <<endl;
    cout<< "Energy variation: "<< S.Energie(N_dim, nu, E) - E0 << endl;
    //Time step
    dt = 2. * pow(10., -7.); //S.pas_temps(t,T,cfl, E, nu, rho); //
    //Computation of forces
    S.Forces(N_dim, nu, E, dt, t , T);
    //Integration of the particle dynamics with the MEMM scheme
    S.Solve(dt,flag_2d, Amortissement, t, T); //Ajouter ici (dans le calcul des vitesses), l'amortissement ?
    //Update of time
    t+= dt;
  }

  //Output of the final solid state
  S.Impression(kimp,rec);

  //Final output
  cout << "Final time of the simulation: " << t<<endl;
  cout << "Computational time: " << (double) (end-start)/CLOCKS_PER_SEC << endl; 
  cout << "Energy variation: " << S.Energie(N_dim, nu, E) - E0 << endl;
  return 0;
}
