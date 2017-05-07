/******************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
******************************************/

template <typename T>
PDEMatrix<T>::PDEMatrix(const long n)
: m_N(n), m_size(pow(n - 1, 2)), m_zero(0)
{}

template <typename T>
PDEMatrix<T>::PDEMatrix(const PDEMatrix<T>& copy)
: m_N(copy.GetN()), m_size(copy.GetSize()), m_zero(copy.GetTolerance())
{}

template <typename T>
PDEMatrix<T>::PDEMatrix(const Matrix<T> * copy)
: m_size(copy->GetSize()), m_zero(copy->GetTolerance())
{
  m_N = sqrt(m_size) + 1;
}

template <typename T>
PDEMatrix<T>::~PDEMatrix()
{} 

template <typename T>
void PDEMatrix<T>::SetTolerance(const T tol)
{
  if (tol < 0)
    throw TolErr("void PDEMatrix<T>::SetTolerance");
  m_zero = tol;
}

template <typename T>
T PDEMatrix<T>::GetTolerance() const
{
  return m_zero;
}

template <typename T>
T PDEMatrix<T>::GetSize() const
{
  return m_size;
}

template <typename T>
long PDEMatrix<T>::GetN() const
{
  return m_N;
}

template <typename T>
Matrix<T>& PDEMatrix<T>::Transpose()
{
  return *this;
}

template <typename T>
bool PDEMatrix<T>::IsDiagDom() const
{
  return true;
}

template <typename T>
T PDEMatrix<T>::operator()(const long row, const long col) const
{
  if (row < 0 || row >= m_size)
    throw SubscriptErr(row, "row Matrix<T>& PDEMatrix<T>::operator()");
  if (col < 0 || col >= m_size)
    throw SubscriptErr(col, "col Matrix<T>& PDEMatrix<T>::operator()");

  if (row > col)
  {
    //calculate lower off diagonal
    if (row == col + 1)
    {
      if (row % (m_N - 1) == 0)
        return 0;
      else
        return -0.25;
    }

    //further diagonals
    if (row == col + m_N - 1)
      return -0.25;
  }
  else
  {
    //calculate upper off diagonal
    if (col == row + 1)
    {
      if (col % (m_N - 1) == 0)
        return 0;
      else
        return -0.25;
    }

    //further diagonals
    if (col == row + m_N - 1)
      return -0.25;

    //calculate diagonal
    if (row == col)
      return 1;
  }

  //everything else
  return 0;
}

template <typename T>
Matrix<T>& PDEMatrix<T>::operator()(const long row __attribute__((unused)), const long col __attribute__((unused)), const T val __attribute__((unused)))
{
  throw PDEErr("Matrix<T>& PDEMatrix<T>::operator()");
}

template <typename T>
unique_ptr<Matrix<T>> PDEMatrix<T>::operator*(const T rhs __attribute__((unused))) const
{
  throw PDEErr("unique_ptr<Matrix<T>> PDEMatrix<T>::operator*");
}

template <typename T>
Vector<T> PDEMatrix<T>::operator*(const Vector<T>& rhs) const
{
  Vector<T> vect = Vector<T>(m_size);
  for (long i = 0; i < m_size; i++)
  {
    T sum = 0;
    for (long j = 0; j < m_size; j++)
      sum += ((*this)(i, j) * rhs[j]);
    vect[i] = sum;
  }
  return vect;
}

template <typename T>
unique_ptr<Matrix<T>> PDEMatrix<T>::operator*(const Matrix<T>& rhs __attribute__((unused))) const
{
  throw PDEErr("unique_ptr<Matrix<T>> PDEMatrix<T>::operator*");
}

template <typename T>
unique_ptr<Matrix<T>> PDEMatrix<T>::operator+(const T rhs __attribute__((unused))) const
{
  throw PDEErr("unique_ptr<Matrix<T>> PDEMatrix<T>::operator+");
}

template <typename T>
unique_ptr<Matrix<T>> PDEMatrix<T>::operator+(const Matrix<T>& rhs __attribute__((unused))) const
{
  throw PDEErr("unique_ptr<Matrix<T>> PDEMatrix<T>::operator+");
}

template <typename T>
Matrix<T>& PDEMatrix<T>::operator+=(const T rhs)
{
  *this = *(*this + rhs);
  return *this;
}

template <typename T>
Matrix<T>& PDEMatrix<T>::operator+=(const Matrix<T>& rhs)
{
  *this = *(*this + rhs);
  return *this;
}

template <typename T>
unique_ptr<Matrix<T>> PDEMatrix<T>::operator-(const T rhs) const
{
  unique_ptr<Matrix<T>> m(new PDEMatrix<T>(*this));
  *m += (rhs * -1);
  return m;
}

template <typename T>
unique_ptr<Matrix<T>> PDEMatrix<T>::operator-(const Matrix<T>& rhs) const
{
  unique_ptr<Matrix<T>> m(new PDEMatrix<T>(*this));
  *m += *(rhs * -1);
  return m;
}

template <typename T>
Matrix<T>& PDEMatrix<T>::operator-=(const T rhs)
{
  *this += (rhs * -1);
  return *this;
}

template <typename T>
Matrix<T>& PDEMatrix<T>::operator-=(const Matrix<T>& rhs)
{
  *this += *(rhs * -1);
  return *this;
}

template <typename T>
Matrix<T>& PDEMatrix<T>::operator=(const Matrix<T>& rhs __attribute__((unused)))
{
  throw PDEErr("Matrix<T>& PDEMatrix<T>::operator=");
}

template <typename T> 
PDEMatrix<T>& PDEMatrix<T>::operator=(const PDEMatrix<T>& rhs)
{
  m_N = rhs.GetN();
  m_size = rhs.GetSize();
  m_zero = rhs.GetTolerance();
  return *this;
}
