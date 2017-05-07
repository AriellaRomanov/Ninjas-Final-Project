/******************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
******************************************/

Jacobi::Jacobi(){}

template <typename T>
Vector<T> Jacobi::operator()(Matrix<T>& matrix, const Vector<T>& b_vect)
{
  if (!matrix.IsDiagDom())
    throw DiagDomErr();
  VNorm1 norm;
  Vector<T> x_vect0 = b_vect;
  Vector<T> x_vect1 = b_vect;
  while (norm((matrix * x_vect0) - b_vect) > matrix.GetTolerance())
  {
    for (long i = 0; i < matrix.GetSize(); i++)
    {
      T sigma = 0;
      for (long j = 0; j < matrix.GetSize(); j++)
      {
        if (i != j)
          sigma += (matrix(i, j) * x_vect0[j]);
      }
      x_vect1[i] = (1 / matrix(i, i)) * (b_vect[i] - sigma);
    }
    x_vect0 = x_vect1;
  }

  return x_vect0;
}

template <typename T>
Vector<T> Jacobi::operator()(PDEMatrix<T>& matrix, const Vector<T>& b_vect)
{
  if (!matrix.IsDiagDom())
    throw DiagDomErr();
  if (matrix.GetSize() != b_vect.GetSize())
    throw SizeErr(matrix.GetSize(), b_vect.GetSize(), "Vector<T> Jacobi::operator()");

  Vector<T> x_vect0(b_vect);
  Vector<T> vect(b_vect);

  Vector<T> err(b_vect);
  for (long i = 0; i < err.GetSize(); i++)
    err[i] = matrix.GetTolerance();

  while (!((matrix * x_vect0) - b_vect < err))
  {
    for (long i = 0; i < matrix.GetSize(); i++)
    {
      T sum = 0;
      if(i + 1 < matrix.GetSize())
        sum += (matrix(i, i + 1) * x_vect0[i + 1]);
      if(i - 1 >= 0)
        sum += (matrix(i, i - 1) * x_vect0[i - 1]);
      if(i - (matrix.GetN() - 1) >= 0)
        sum += (matrix(i, i - (matrix.GetN() - 1)) * x_vect0[i - (matrix.GetN() - 1)]);
      if(i + (matrix.GetN() - 1) < matrix.GetSize())
        sum += (matrix(i, i + (matrix.GetN() - 1)) * x_vect0[i + (matrix.GetN() - 1)]);
      vect[i] = sum;
    }
    x_vect0 = b_vect - vect;
  }

  return x_vect0;
}