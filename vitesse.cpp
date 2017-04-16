#include "geometry.hpp"

//Initial velocity of the solid particles
Vector_3 velocity(const Point_3 &p)
{
  return Vector_3(0.01*(p.x()-1.),0,0);
  //return Vector_3(0,0,0);
}

//Initial angular velocity of the solid particles
Vector_3 omega(const Point_3 &p)
{
  return Vector_3(0,0,0);
  //return Vector_3(0.1,0,0);
  
}
