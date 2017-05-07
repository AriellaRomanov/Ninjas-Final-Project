/******************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
******************************************/

template <typename T>
PDESolution<T>::PDESolution(T x0(T), T x1(T), T y0(T), T y1(T), T g(T))
{
  m_x0 = x0;
  m_x1 = x1;
  m_y0 = y0;
  m_y1 = y1;
  m_g = g; 
}

template <typename T>
Vector<T> PDESolution<T>::Gaussian(const long N, const T tolerance) const
{
  long size = N - 1;
  Vector<T> b_vector(size * size);
  b_vector.SetTolerance(tolerance);
  for(int i = 0; i < size * size; i++)
    b_vector[i] = 0;

  for(int i = 0; i < size; i++)
  {
    b_vector[i] += m_x0(static_cast<T>(i + 1) / N) / 4;
    b_vector[i * size] += m_y0(static_cast<T>(i + 1) / N) / 4;
    b_vector[((size * size - size) + i)] += m_x1(static_cast<T>(i + 1) / N) / 4;
    b_vector[((i + 1) * size) - 1] += m_y1(static_cast<T>(i + 1) / N) / 4;
    for(int k = 0; k < size; k++)
      b_vector[(i * size) + k] -= m_g(static_cast<T>(i + 1) / N) * (1 / static_cast<T>(N * N * 4));
  }

  class Gaussian gauss;
  DenseMatrix<T> A(size * size);
  A.SetTolerance(tolerance);

  for(long row = 0; row < size * size; row++)
  {
    for(long col = 0; col < size * size; col++)
    {
      A(row, col, 0);
      //calculate diagonal
      if (row == col)
        A(row, col, 1);

      //calculate lower off diagonal
      else if (row == col + 1 && row % size != 0)
        A(row, col, -0.25);

      //calculate upper off diagonal
      else if (col == row + 1 && col % size != 0)
        A(row, col, -0.25);

      //further diagonals
      else if (row == col + size || col == row + size)
        A(row, col, -0.25);   
    }
  }

  return gauss(A, b_vector);
}

template <typename T>
Vector<T> PDESolution<T>::Jacobi(const long N, const T tolerance) const
{
  long size = N - 1;
  Vector<T> b_vector(size * size);
  for(int i = 0; i < (size * size); i++)
    b_vector[i] = 0;
  for(int i = 0; i < size; i++)
  {
    b_vector[i] += m_x0(static_cast<T>(i + 1) / N) / 4;
    b_vector[i * size] += m_y0(static_cast<T>(i + 1) / N) / 4;
    b_vector[(size * size) - size + i] += m_x1(static_cast<T>(i + 1) / N) / 4;
    b_vector[((i + 1) * size) - 1] += m_y1(static_cast<T>(i + 1) / N) / 4;
    
    for(int k = 0; k < size; k++)
      b_vector[(i * size) + k] -= m_g(static_cast<T>(i + 1) / N) * (1 / static_cast<T>(N * N * 4));
  }

  PDEMatrix<T> a_matrix(N);
  a_matrix.SetTolerance(tolerance);
  b_vector.SetTolerance(tolerance);
  
  class Jacobi jacob;
  return jacob(a_matrix, b_vector);
}