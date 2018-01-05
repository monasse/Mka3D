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
\authors Laurent Monasse
 *  \file forces_ext.cpp
 *  \brief Definition of the functions giving the external forces.
 */

#include "geometry.hpp"
#include "solide.hpp"
#ifndef FORCES_EXT_CPP
#define FORCES_EXT_CPP

//Pression dans tube
Vector_3 Forces_externes(const double& t, const double& T, const Face& face, const double& mu, const int& fixe)
{
  double p_max = 5. * 90000000. * 0.75 * 2. * mu; //En Pa pression max, 5 fois seuil elas théorique...
  //double p = p_max * t / T;
  double p = p_max / 50.;

  Point_3 pos_centre = face.centre;
  double r = sqrt( pos_centre.x() * pos_centre.x() + pos_centre.y() * pos_centre.y());

  //Changer la face de reconnaitre les faces !!!
  //Faire avec des BC de type fixe == 2 ou 3 et remplir ça en vérifiant que particules ont des voisins ou pas !
  if(fixe == 2)
    return p * (-face.normale);
  else
    return Vector_3(0., 0., 0.);
}

#endif
