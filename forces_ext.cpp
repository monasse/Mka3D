/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 

Contributors: Laurent Monasse
*/

/*!
 *  \file forces_ext.cpp
 *  \brief Definition of the functions giving the external forces.
 */

#include "geometry.hpp"
#ifndef FORCES_EXT_CPP
#define FORCES_EXT_CPP

Vector_3 Forces_externes(const Point_3 &X, const Vector_3 &e)
{
  return Vector_3(0,0,0);
}

Vector_3 Moments_externes(const Point_3 &X, const Vector_3 &e)
{
  return Vector_3(0,0,0);
}


#endif
