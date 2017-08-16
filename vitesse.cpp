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

#include "geometry.hpp"

//Initial velocity of the solid particles
Vector_3 velocity(const Point_3 &p)
{
  //return Vector_3(0.01*(p.x()-1.),0,0);
  if(p.x() <= 1.)
     return Vector_3(-0.01,0,0);
  else if(p.x() >= 19.)
    return Vector_3(0,0,0);
}

//Initial angular velocity of the solid particles
Vector_3 omega(const Point_3 &p)
{
  return Vector_3(0,0,0);
  //return Vector_3(0.1,0,0);
  
}
