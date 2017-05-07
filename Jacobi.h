/******************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
******************************************/

#pragma once

#include "Matrix.h"
#include "VNorm1.h"

class Jacobi
{
  public:
    // Desc: Constructor for Jacobi class
    // Pre: none
    // Post: instantiates an instance of the Jacobi class
    Jacobi();
    // Desc: applies the Jacobi process to a Matrix of type T
    // Pre: 
    // Post: 
    template <typename T>
    Vector<T> operator()(Matrix<T>& matrix, const Vector<T>& b_vect);
    // Desc: applies the Jacobi process to a PDEMatrix of type T
    // Pre: 
    // Post: 
    template <typename T>
    Vector<T> operator()(PDEMatrix<T>& matrix, const Vector<T>& b_vect);
};

#include "Jacobi.hpp"