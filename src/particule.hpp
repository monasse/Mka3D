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

#ifndef PARTICULE_HPP
#define PARTICULE_HPP

#include "geometry.hpp"
#include "solide.hpp"
#include "face.hpp"

class Solide;

class Face;

class Particule
{

 public:
  //Functions
  Particule();//:faces(std::vector<Face>(1)){}
  Particule(const int& Id); // Permet de créer une particule à partir de son numéro et de la donnée de ses vertex
  Particule(const double &x_min, const double &y_min, const double &z_min, 
			const double &x_max, const double &y_max,const double &z_max);
  Particule(const Point_3 &c, const double &x_min, const double &y_min, const double &z_min, 
			const double &x_max, const double &y_max,const double &z_max, 
			const std::vector<Face> & F);
  ~Particule();
  Particule & operator=(const Particule &P); // opérateur = surcharge pour l'affectation
  void solve_position(const double &dt, const bool &flag_2d, const double& t, const double& T);
  void solve_vitesse(const double &dt, const bool &flag_2d, const double& Amort, const double& t, const double& T);

  Vector_3 vitesse_parois(const Point_3& X_f);  
  Vector_3 vitesse_parois_prev(const Point_3& X_f);

  void barycentre(Solide* Sol, const int& cell_type); //Calcul le barycentre d'une particule
  void volume(Solide* Sol, const int& cell_type); //calcul le volume d'une particule
  bool contient_face(const Face& f); //Renvoie vraie si particule contient les 3 vertex de la face
  void calcul_diametre(Solide* Sol);

  //Attributs
  std::vector<int> faces; //!< liste de faces de la particule
  std::vector<int> vertices; //Utile pour retrouver les faces lors de l'importation de la connectivité
  int fixe; //!< =true si la particule est fixe, false sinon
  double m; //!< Masse de la particule
  double V; //!< Volume de la particule
  Point_3 x0; //!<Position du centre de la particule dans la configuration initiale
  Vector_3 Dx; //!<Depleacement du centre de la particule en t
  //Vector_3 Dx_plas; //Déplacements plastiques...
  Vector_3 Dxprev; //!<Deplacement du centre de la particule en t-dt
  Vector_3 Fi; //!<Forces int&eacute;rieures du solide
  Vector_3 u; //!< Vitesse de la particule au temps t
  Vector_3 u_prev; //!< Vitesse de la particule au temps t-dt/2
  Aff_transformation_3 mvt_t; //!<Transformation affine de la particule au temps t
  Aff_transformation_3 mvt_tprev; //!<Transformation affine de la particule au temps t-dt
  double h; //Diamètre de la particule

  //Variables pour plasticité et nouvelle formulation Mka !
  Matrix grad; //Gradient reconstruit par particule
  Matrix discrete_gradient; //Gradient symétrique reconstruit par particule
  Matrix contrainte; //Contrainte par particule
  double def_plas_cumulee; //Déformation plastique cumulée du lien
  Matrix epsilon_p; //Déformation plastique rémanante
  double seuil_elas;

  //Pour compensation des sommations
  Vector_3 err_u;
  Vector_3 err_Dx;

  int id; //Numéro de la particule dans la map de Solide
  int BC; //Vaut 0 si particule pas au bord et peut valoir 1,2,etc... selon la condition de bord

  //Pour indiquer s'il faut prendre en compte ou pas la particule car elle a été splittée
  bool split;
}; 

#endif
