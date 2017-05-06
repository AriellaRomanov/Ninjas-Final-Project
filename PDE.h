#pragma once
/* #include "PDEMatrix.h"
#include "DenseMatrix.h" */

/* template <typename T, T x0(T), T x1(T), T y0(T), T y1(T), T g(T)>
Vector<T> PDESolution(long method); */


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
    PDESolution(T x0(T), T x1(T), T y0(T), T y1(T), T g(T));
    Vector<T> Gaussian(const long N, const T tolerance) const;
    Vector<T> Jacobi(const long N, const T tolerance) const;
};

#include "PDE.hpp"