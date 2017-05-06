#pragma once
/* #include "PDEMatrix.h"
#include "DenseMatrix.h" */

/* template <typename T, T x0(T), T x1(T), T y0(T), T y1(T), T g(T)>
Vector<T> PDESolution(long method); */




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
    PDESolution(T x0(T), T x1(T), T y0(T), T y1(T), T g(T));
    Vector<T> Gaussian(long m_n);
    Vector<T> Jacobi(long iter_max, long m_n);
};

#include "PDE.hpp"