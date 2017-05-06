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
    // Pre: 
    // Post: 
    template <typename T>
    Vector<T> operator()(Matrix<T>& matrix, const Vector<T>& b_vect);
};

#include "Gaussian.hpp"