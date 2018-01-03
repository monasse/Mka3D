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
 *  \file solide.hpp
 *  \brief D&eacute;finition des classes d&eacute;crivant le solide.
 */

#ifndef SOLIDE_HPP
#define SOLIDE_HPP

#include <map>
#include "geometry.hpp"
#include "particule.hpp"

//! Definition de la classe Solide  
class Solide
{
	
public:
  
  Solide();//:solide(std::vector<Particule>(1)){}
  Solide(const double& E, const double& nu);
  Solide(const std::vector<Particule> & Part);
  ~Solide();
  Solide & operator=(const Solide &S); // opérateur = surcharge pour l'affectation
  //void Affiche();  //fonction auxilaire utile pour les test
  int size(){
	return solide.size();
  }
  void Impression(const int &n, const bool &reconstruction);
  void Init(const char* s1, const char* s2, const char* s3, const bool& rep, const int& numrep, const double& rho);
  void Solve_position(const double &dt, const bool &flag_2d, const double& t, const double& T);
  //void stock_def_plastique(const double &dt);
  void Solve_vitesse(const double &dt, const bool &flag_2d, const double& Amort, const double& t, const double& T);
  void Forces(const int &N_dim, const double &nu, const double &E, const double& dt, const double& t, const double& T);
  void Forces_internes(const int &N_dim, const double &nu, const double &E, const double& dt);
  void update_triangles();
  //void breaking_criterion();
  double Energie(const int &N_dim, const double &nu, const double &E);
  double Energie_potentielle(const int &N_dim, const double &nu, const double &E);
  double Energie_cinetique();
  double pas_temps(const double &t, const double &T, const double &cfls, const double &E, const double &nu, const double &rhos);
  // private :
  //std::vector<Particule> solide; //!< Maillage solide
  std::map<int, Particule> solide;

  double lambda; //Premier coeff de lamé
  double mu; //Second coefficient de lamé
};


#endif
