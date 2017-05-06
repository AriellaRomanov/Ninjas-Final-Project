#pragma once

#include <cmath>
#include "Vector.h"

class VNorm1
{
  public:
    // Desc: Constructor for VNorm1 class
    // Pre: None
    // Post: Instantiates an instance of the VNorm1 class
    VNorm1();
    // Desc: applies the 1 norm to a Vector of type T
    // Pre: Type T must have defined for it the absolute value
    // operator and the addition operator
    // Post: Returns the sum of the absolute value of the elements of the Vector
    template <typename T>
    T operator()(const Vector<T>& vect);
};

#include "VNorm1.hpp"