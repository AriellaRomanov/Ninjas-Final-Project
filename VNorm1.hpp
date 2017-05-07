/******************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
******************************************/

VNorm1::VNorm1() {}

template <typename T>
T VNorm1::operator()(const Vector<T>& vect)
{
  T result = 0;
  for (long i = 0; i < vect.GetSize(); i++)
    result += abs(vect[i]);
  return result;
}