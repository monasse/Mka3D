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
 *  \file geometry.cpp
 *  \brief Definition of the classes involved in geometry computation.
 */

#include <iostream>
#include "geometry.hpp"

#ifndef GEOMETRY_CPP
#define GEOMETRY_CPP

//////////////////////////////////////////////
//  Class Bbox                              //
//////////////////////////////////////////////

//Default constructor
Bbox::Bbox()
{
  xm = ym = zm = std::numeric_limits<double>::infinity();
  xM = yM = zM = -std::numeric_limits<double>::infinity();
}

Bbox::Bbox(const double& x_min, const double& y_min, const double& z_min, const double& x_max, const double& y_max, const double& z_max)
{
  xm = x_min;
  ym = y_min;
  zm = z_min;
  xM = x_max;
  yM = y_max;
  zM = z_max;
}

double Bbox::xmin() const
{
  return xm;
}

double Bbox::ymin() const
{
  return ym;
}

double Bbox::zmin() const
{
  return zm;
}

double Bbox::xmax() const
{
  return xM;
}

double Bbox::ymax() const
{
  return yM;
}

double Bbox::zmax() const
{
  return zM;
}

Bbox Bbox::operator+(const Bbox &bb) const
{
  return Bbox(min(xm,bb.xmin()),min(ym,bb.ymin()),min(zm,bb.zmin()),max(xM,bb.xmax()),max(yM,bb.ymax()),max(zM,bb.zmax()));
}

Bbox& Bbox::operator+=(const Bbox &bb)
{
  xm = min(xm,bb.xmin());
  ym = min(ym,bb.ymin());
  zm = min(zm,bb.zmin());
  xM = max(xM,bb.xmax());
  yM = max(yM,bb.ymax());
  zM = max(zM,bb.zmax());
}

bool do_overlap(const Bbox &bb1, const Bbox &bb2)
{
  if (bb1.xmax() < bb2.xmin() || bb2.xmax() < bb1.xmin())
    return false;
  if (bb1.ymax() < bb2.ymin() || bb2.ymax() < bb1.ymin())
    return false;
  
  return true;
}



//////////////////////////////////////////////
/// Class Point_3                         ////
//////////////////////////////////////////////

//Default constructor
Point_3::Point_3()
{
  p[0] = p[1] = p[2] = 0.;
}

//Constructor from coordinates
Point_3::Point_3(const double& x0, const double& y0, const double& z0)
{
  p[0] = x0;
  p[1] = y0;
  p[2] = z0;
}

//Extraction of index i
double Point_3::operator[](const int& i) const
{
  return p[i];
}

double Point_3::x() const
{
  return p[0];
}

double Point_3::y() const 
{
  return p[1];
}

double Point_3::z() const
{
  return p[2];
}


Bbox Point_3::bbox() const
{
  return Bbox(p[0],p[1],p[2],p[0],p[1],p[2]);
}

ostream& operator<<(ostream &os, const Point_3 &p)
{
  os << p.x() << " " << p.y() << " " << p.z();
  return os;
}


double squared_distance(const Point_3 &p, const Point_3 &q)
{
  return Vector_3(p,q).squared_length();
}



//////////////////////////////////////////
//    Class Vector_3                   ///
//////////////////////////////////////////

//Default constructor
Vector_3::Vector_3()
{
  vec[0] = vec[1] = vec[2] = 0;
}

//Constructor from coordinates
Vector_3::Vector_3(const double& x0, const double& y0, const double& z0)
{
  vec[0] = x0;
  vec[1] = y0;
  vec[2] = z0;
}

//Constructor from points
Vector_3::Vector_3(const Point_3& p1, const Point_3& p2)
{
  vec[0] = p2[0]-p1[0];
  vec[1] = p2[1]-p1[1];
  vec[2] = p2[2]-p1[2];
}

//Extraction of index i
double Vector_3::operator[](const int& i) const
{
  return vec[i];
}

double Vector_3::x() const
{
  return vec[0];
}

double Vector_3::y() const 
{
  return vec[1];
}

double Vector_3::z() const
{
  return vec[2];
}

//Opposite vector
Vector_3 Vector_3::operator-() const
{
  return Vector_3(-vec[0],-vec[1],-vec[2]);
}

//Division by a scalar
Vector_3 Vector_3::operator/(const double &s) const
{
  if(s!=0){
    return Vector_3(vec[0]/s,vec[1]/s,vec[2]/s);
  }
  else {
    throw std::invalid_argument( "Division by zero" );
  }
}

//Squared length of vector
double Vector_3::squared_length() const
{
  return vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
}

//Scalar product
double Vector_3::operator*(const Vector_3 &v) const
{
  return vec[0]*v[0]+vec[1]*v[1]+vec[2]*v[2];
}

//Right multiplication by a scalar
Vector_3 Vector_3::operator*(const double &s) const
{
  return Vector_3(s*vec[0],s*vec[1],s*vec[2]);
}

//Left multiplication by a scalar
Vector_3 operator*(const double &s, const Vector_3 &v)
{
  return Vector_3(s*v[0],s*v[1],s*v[2]);
}

ostream& operator<<(ostream &os, const Vector_3 &v)
{
  os << v.x() << " " << v.y() << " " << v.z();
  return os;
}

//Addition of vectors
Vector_3 Vector_3::operator+(const Vector_3 &v2) const
{
  return Vector_3(vec[0]+v2[0],vec[1]+v2[1],vec[2]+v2[2]);
}

//Subtraction of vectors
Vector_3 Vector_3::operator-(const Vector_3 &v2) const
{
  return Vector_3(vec[0]-v2[0],vec[1]-v2[1],vec[2]-v2[2]);
}

//Addition of a vector to a point
Point_3 operator+(const Point_3 &p, const Vector_3 &v)
{
  return Point_3(p[0]+v[0],p[1]+v[1],p[2]+v[2]);
}

//Subtraction of a vector to a point
Point_3 operator-(const Point_3 &p, const Vector_3 &v)
{
  return Point_3(p[0]-v[0],p[1]-v[1],p[2]-v[2]);
}

//Vector difference of two points
Vector_3 operator-(const Point_3 &p1, const Point_3 &p2)
{
  return Vector_3(p1,p2);
}

//////////////////////////////////////////////////////////
//  Affine transformations class                        //
//////////////////////////////////////////////////////////

//Default constructor: identity
Aff_transformation_3::Aff_transformation_3()
{
  lin[0][0] = lin[1][1] = lin[2][2] = 1;
  lin[0][1] = lin[0][2] = lin[1][0] = lin[1][2] = lin[2][0] = lin[2][1] = 0;
  translation = Vector_3(0,0,0);
}

//Constructor from the matrix coefficients m
Aff_transformation_3::Aff_transformation_3(const double& m00, const double& m10, const double& m20, const double& m01, const double& m11, const double& m21, const double& m02, const double& m12, const double& m22, const double& m03, const double& m13, const double& m23)
{
  lin[0][0] = m00;
  lin[1][0] = m10;
  lin[2][0] = m20;
  lin[0][1] = m01;
  lin[1][1] = m11;
  lin[2][1] = m21;
  lin[0][2] = m02;
  lin[1][2] = m12;
  lin[2][2] = m22;
  translation = Vector_3(m03,m13,m23);
}

//Constructor for a pure translation
Aff_transformation_3::Aff_transformation_3(const Vector_3 &v)
{
  lin[0][0] = lin[1][1] = lin[2][2] = 1.;
  lin[1][0] = lin[2][0] = lin[0][1] = lin[2][1] = lin[0][2] = lin[1][2] = 0.;
  translation = v;
}

//Constructor from the matrix coefficients m without translation
Aff_transformation_3::Aff_transformation_3(const double& m00, const double& m10, const double& m20, const double& m01, const double& m11, const double& m21, const double& m02, const double& m12, const double& m22)
{
  lin[0][0] = m00;
  lin[1][0] = m10;
  lin[2][0] = m20;
  lin[0][1] = m01;
  lin[1][1] = m11;
  lin[2][1] = m21;
  lin[0][2] = m02;
  lin[1][2] = m12;
  lin[2][2] = m22;
  translation = Vector_3(0,0,0);
}

Aff_transformation_3 Aff_transformation_3::operator*(const Aff_transformation_3 &T) const
{
  double Tp[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Tp[i][j] = 0.;
      for(int k=0;k<3;k++){
	Tp[i][j] += lin[i][k]*T.lin[k][j];
      }
    }
  }
  Vector_3 v = transform(T.translation)+translation;
  return Aff_transformation_3(Tp[0][0],Tp[1][0],Tp[2][0],Tp[0][1],Tp[1][1],Tp[2][1],Tp[0][2],Tp[1][2],Tp[2][2],v[0],v[1],v[2]);
}


//Transformation of a point
Point_3 Aff_transformation_3::transform(const Point_3 &p) const
{
  double coord[3];
  for(int i=0;i<3;i++){
    coord[i] = 0.;
    for(int j=0;j<3;j++){
      coord[i] += lin[i][j]*p[j];
    }
  }
  return Point_3(coord[0],coord[1],coord[2])+translation;
}

Point_3 Aff_transformation_3::operator()(const Point_3 &p) const
{
  return transform(p);
}


//Transformation of a vector
Vector_3 Aff_transformation_3::transform(const Vector_3 &v) const
{
  double coord[3];
  for(int i=0;i<3;i++){
    coord[i] = 0.;
    for(int j=0;j<3;j++){
      coord[i] += lin[i][j]*v[j];
    }
  }
  return Vector_3(coord[0],coord[1],coord[2]);
}

Vector_3 Aff_transformation_3::operator()(const Vector_3 &v) const
{
  return transform(v);
}





///////////////////////////////////////////////////////////
///  Miscellanious geometric functions                  ///
///////////////////////////////////////////////////////////


//Cross product
Vector_3 cross_product(const Vector_3 &v1, const Vector_3 &v2)
{
  return Vector_3(v1[1]*v2[2]-v1[2]*v2[1],v1[2]*v2[0]-v1[0]*v2[2],v1[0]*v2[1]-v1[1]*v2[0]);
}

//Centroid of points
Point_3 centroid(const std::vector<Point_3>::iterator &begin, const std::vector<Point_3>::iterator &end)
{
  double x,y,z;
  x = y = z = 0;
  int n=0;
  for(std::vector<Point_3>::iterator it=begin;it!=end;++it){
    x += (*it)[0];
    y += (*it)[1];
    z += (*it)[2];
    n += 1;
  }
  if(n>0){
    return Point_3(x/n,y/n,z/n);
  } else {
    throw std::invalid_argument( "Division by zero" );
  }
}

//Vector normal to the plane defined by three points (in direct order with regards to the normal)
Vector_3 orthogonal_vector(const Point_3 &p0, const Point_3 &p1, const Point_3 &p2)
{
  Vector_3 v1(p0,p1);
  Vector_3 v2(p0,p2);
  Vector_3 n = cross_product(v1,v2);
  if(n.squared_length()==0){
    cout << "p0=(" << p0 << "), p1=(" << p1 << "), p2=(" << p2 << ")"<< endl;
    cout << "p1-p0=(" << Vector_3(p0,p1) << ")" << endl;
    cout << "p2-p0=(" << Vector_3(p0,p2) << ")" << endl;
    cout << "n=(" << n << ")"<< endl;
    throw std::invalid_argument( "Colinear points" );
  } else {
    return n;
  }
}

//////////////////////////////////////////////
// class Triangle_3                         //
//////////////////////////////////////////////

//Default constructor
Triangle_3::Triangle_3()
{
}

//Constructor overload
Triangle_3::Triangle_3(const Point_3 &p, const Point_3 &q, const Point_3 &r)
{
  v[0] = p;
  v[1] = q;
  v[2] = r;
}

//Vertex number i
Point_3 Triangle_3::vertex(const int& i) const
{
  return v[i];
}

Point_3 Triangle_3::operator[](const int& i) const
{
  return v[i];
}

bool Triangle_3::is_degenerate() const
{
  return cross_product(Vector_3(v[0],v[1]),Vector_3(v[0],v[2])).squared_length()==0.;
}

double Triangle_3::squared_area() const
{
  return cross_product(Vector_3(v[0],v[1]),Vector_3(v[0],v[2])).squared_length()/4.;
}

Bbox Triangle_3::bbox() const
{
  return v[0].bbox()+v[1].bbox()+v[2].bbox();
}

Triangle_3 Triangle_3::transform(const Aff_transformation_3 &T) const
{
  return Triangle_3(T(v[0]),T(v[1]), T(v[2]));
}

///////////////////////////////////////////////////
//    class Tetrahedron                          //
///////////////////////////////////////////////////

//Default constructor
Tetrahedron::Tetrahedron()
{
}

Tetrahedron::Tetrahedron(const Point_3 &p0, const Point_3 &p1, const Point_3 &p2, const Point_3 &p3)
{
  v[0] = p0;
  v[1] = p1;
  v[2] = p2;
  v[3] = p3;
}

Point_3 Tetrahedron::vertex(const int& i) const
{
  return v[i];
}

Point_3 Tetrahedron::operator[](const int& i) const
{
  return v[i];
}

bool Tetrahedron::is_degenerate() const
{
  return cross_product(Vector_3(v[0],v[1]),Vector_3(v[0],v[2]))*Vector_3(v[0],v[3])==0.;
}

double Tetrahedron::volume() const
{
  return cross_product(Vector_3(v[0],v[1]),Vector_3(v[0],v[2]))*Vector_3(v[0],v[3])/6.;
}

Bbox Tetrahedron::bbox() const
{
  return v[0].bbox()+v[1].bbox()+v[2].bbox()+v[3].bbox();
}

Tetrahedron Tetrahedron::transform(const Aff_transformation_3 &T) const
{
  return Tetrahedron(T(v[0]),T(v[1]),T(v[2]),T(v[3]));
}




#endif
