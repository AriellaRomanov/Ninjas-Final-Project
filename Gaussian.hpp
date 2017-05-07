/******************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
******************************************/

Gaussian::Gaussian(){}

template <typename T>
Vector<T> Gaussian::operator()(Matrix<T>& matrix, const Vector<T>& b_vect)
{
  if (matrix.GetSize() != b_vect.GetSize())
    throw SizeErr(matrix.GetSize(), b_vect.GetSize(), "Vector<T> Gaussian::operator()");

  long size = matrix.GetSize();

  for(long i = 0; i < size; i++)
  {
    for(long j = i + 1; j < size;j++)
    {
      if(matrix(i, j) != 0)
      {
        T temp = matrix(i, j) / matrix(i, i);
        for(long k = i; k < size; k++)
          matrix(k, j, matrix(k, j) - (temp * matrix(k, i)));
        b_vect[j] -= temp * b_vect[i];
      }
    }
  }

  //Bottom half is 0's
  for(long i = size - 1; i >= 0; i--)
  {
    for(long k = 0; k < i; k++)
    {
      b_vect[k] -= (matrix(i, k) / matrix(i, i)) * b_vect[i];
      matrix(i, k, 0);
    }
  }
  
  //Top half is 0's
  for(long i = 0; i < size; i++)
  {
    if (matrix(i, i) == 0)
      throw  DivByZeroErr("Gaussian");
    b_vect[i] = b_vect[i] / matrix(i, i);
  }

  //Diagonal is 1's
  return b_vect;
}