/******************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
******************************************/

#pragma once

#include "Jacobi.h"
#include "Gaussian.h"

template <typename T>
class PDESolution
{
  private:
    T (*m_x0)(T);
    T (*m_x1)(T);
    T (*m_y0)(T);
    T (*m_y1)(T);
    T (*m_g)(T);
  public:
    // Desc: Instantiates a new matrix of the given size
    // Pre: Size must be > 0
    // Post: A new matrix will be created of the given dimension
    PDESolution(T x0(T), T x1(T), T y0(T), T y1(T), T g(T));
    // Desc: Creates a PDEMatrix and a vector to solve Ax=b
    // using Gaussian elimination
    // Pre: N must be greater than one, tolerance
    // must be greater than zero
    // Post: A vector solution of the given callbacks will
    // be returned
    Vector<T> Gaussian(const long N, const T tolerance) const;
    // Desc: Creates a PDEMatrix and a vector to solve Ax=b
    // using Jacobi elimination
    // Pre: N must be greater than one, tolerance
    // must be greater than zero
    // Post: A vector solution of the given callbacks will
    // be returned
    Vector<T> Jacobi(const long N, const T tolerance) const;
};

#include "PDE.hpp"