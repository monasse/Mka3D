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
 *  \file geometry.hpp
 *  \brief D&eacute;finition des classes d&eacute;crivant la g&eacute;om&eacute;trie.
 */

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <vector>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <cassert>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdlib.h>

using namespace std;


class Bbox
{
private:
  double xm,ym,zm,xM,yM,zM;
public:
  Bbox();
  Bbox(const double &x_min, const double &y_min, const double &z_min, const double &x_max, const double &y_max, const double &z_max);
  double xmin() const;
  double ymin() const;
  double zmin() const;
  double xmax() const;
  double ymax() const;
  double zmax() const;
  Bbox operator+(const Bbox &bb) const;
  Bbox& operator+=(const Bbox &bb);
};

bool do_overlap(const Bbox &bb1, const Bbox &bb2);


class Point_3 
{
private:
  double p[3];
public:
  Point_3();
  Point_3(const double &x0, const double &y0, const double &z0);
  double operator[](const int &i) const;
  Point_3 operator/(const double &s) const;
  double x() const;
  double y() const;
  double z() const;
  Bbox bbox() const;
};

ostream& operator<<(ostream &os, const Point_3 &p);
double squared_distance(const Point_3 &p, const Point_3 &q);
bool operator==(const Point_3 &vec1, const Point_3 &vec2);
Point_3 operator+(const Point_3 &p, const Point_3 &v);

class Vector_3 
{
private :
  double vec[3];
public:
  Vector_3();
  Vector_3(const double &x0, const double &y0, const double &z0);
  Vector_3(const Point_3 &p1, const Point_3 &p2);
  double operator[](const int &i) const;
  double x() const;
  double y() const;
  double z() const;
  Vector_3 operator-() const;
  Vector_3 operator/(const double &s) const;
  double squared_length() const;
  double operator*(const Vector_3 &v) const;
  Vector_3 operator*(const double &s) const;
  Vector_3 operator+(const Vector_3 &v2) const;
  Vector_3 operator-(const Vector_3 &v2) const;
};

ostream& operator<<(ostream &os, const Vector_3 &v);
Point_3 operator+(const Point_3 &p, const Vector_3 &v);
Point_3 operator-(const Point_3 &p, const Vector_3 &v);
Vector_3 operator-(const Point_3 &p1, const Point_3 &p2);
Vector_3 operator*(const double &s, const Vector_3 &v);
bool operator==(const Vector_3 &vec1, const Vector_3 &vec2);

class Matrix {
 public:
  Matrix();
  Matrix(const Vector_3& colonne_1, const Vector_3& colonne_2, const Vector_3& colonne_3);
  //double norme() const;
  double tr(); //Trace d'une matrice
  Vector_3 c1() const;
  Vector_3 c2() const;
  Vector_3 c3() const;
  Matrix T() const; //Transposée
  //void empty(); //Remet tous les coefficients de la matrice à 0.
  Matrix dev(); //Renvoie le deviateur du tenseur considéré
  double VM(); //Renvoie la norme de Von Mises associée à une matrice
  Matrix operator/(const double& rel);
  double norme() const; //Norme 2 au sens des matrices
  Matrix& operator+=(const Matrix &mat);

 public:
  Vector_3 col1;
  Vector_3 col2;
  Vector_3 col3;
};

bool operator==(const Matrix& vec1, const Matrix& vec2);
Matrix operator+(const Matrix& vec1, const Matrix& vec2);
Matrix operator-(const Matrix& vec1, Matrix const& vec2);
Matrix operator-(const Matrix& vec);
Matrix operator*(const double& rel, const Matrix& vec); //Produit scalaire matrice
Matrix operator*(const Matrix& vec, const double& rel); //Produit scalaire matrice
Matrix operator*(const Matrix& vec1, const Matrix& vec2); //Produit simplement contracté
Vector_3 operator*(const Matrix& vec1, const Vector_3& vec2); //Produit matrice vecteur
double contraction_double(const Matrix& vec1, const Matrix & vec2); //Produit doublement contracté
Matrix tens(const Vector_3& vec1, const Vector_3& vec2); //Produit tensoriel
Matrix tens_sym(const Vector_3& vec1, const Vector_3& vec2); //Produit tensoriel symétrique
Matrix unit(); //Matrice unité

class Aff_transformation_3
{
private :
  double lin[3][3];
  Vector_3 translation;
public:
  Aff_transformation_3();
  Aff_transformation_3(const double &m00, const double &m10, const double& m20, const double& m01, const double& m11, const double& m21, const double& m02, const double& m12, const double& m22, const double& m03, const double& m13, const double& m23);
  Aff_transformation_3(const Vector_3 &v);
  Aff_transformation_3(const double& m00, const double& m10, const double& m20, const double& m01, const double& m11, const double& m21, const double& m02, const double& m12, const double& m22);
  Aff_transformation_3 operator*(const Aff_transformation_3 &T) const;
  Point_3 transform(const Point_3 &p) const;
  Vector_3 transform(const Vector_3 &v) const;
  Point_3 operator()(const Point_3 &p) const;
  Vector_3 operator()(const Vector_3 &v) const;
};

  

Vector_3 cross_product(const Vector_3 &v1, const Vector_3 &v2);

Point_3 centroid(const std::vector<Point_3>::iterator &begin, const std::vector<Point_3>::iterator &end);

Vector_3 orthogonal_vector(const Point_3 &p0, const Point_3 &p1, const Point_3 &p2);

class Triangle_3
{
private:
  Point_3 v[3];
public:
  Triangle_3();
  Triangle_3(const Point_3 &p, const Point_3 &q, const Point_3 &r);
  Point_3 vertex(const int& i) const;
  Point_3 operator[](const int& i) const;
  bool is_degenerate() const;
  double squared_area() const;
  Bbox bbox() const;
  Triangle_3 transform(const Aff_transformation_3 &T) const;
};

typedef std::vector<Triangle_3>  Triangles;

class Tetrahedron
{
private:
  Point_3 v[4];
public:
  Tetrahedron();
  Tetrahedron(const Point_3 &p0, const Point_3 &p1, const Point_3 &p2, const Point_3 &p3);
  Point_3 vertex(const int &i) const;
  Point_3 operator[](const int &i) const;
  bool is_degenerate() const;
  double volume() const;
  Bbox bbox() const;
  Tetrahedron transform(const Aff_transformation_3 &T) const;
};

//Function absolute value
double abs(const double &x);



// bool inside_box(const Bbox& cell, const Point_3& P);
// bool box_inside_convex_polygon(const Particule& S, const Bbox& cell);  
// bool inside_convex_polygon(const Particule& S, const Point_3& P);  
// double Error(Solide& S1, Solide& S2);
// void Copy_f_m(Solide& S1, Solide& S2);
// bool box_inside_tetra(const Tetrahedron &tetra, const Bbox& cell);
// bool inside_tetra(const Tetrahedron &tetra, const Point_3& P);

#endif
