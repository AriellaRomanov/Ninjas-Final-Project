/******************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
******************************************/

#pragma once

#include "Matrix.h"
#include "VNorm1.h"

class Gaussian
{
  public:
    // Desc: Constructor for Gaussian class
    // Pre: none
    // Post: instantiates an instance of the Gaussian class
    Gaussian();
    // Desc: applies the Gaussian process to a Matrix of type T
    // Pre: matrix and b_vect must be of the same size
    // Post: A vector solution to the Ax = b problem will be returned
    template <typename T>
    Vector<T> operator()(Matrix<T>& matrix, const Vector<T>& b_vect);
};

#include "Gaussian.hpp"