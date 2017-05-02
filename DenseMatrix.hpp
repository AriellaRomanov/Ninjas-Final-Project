template <typename T>
DenseMatrix<T>::DenseMatrix(const long size)
: m_size(size), m_data(new Vector<T>[m_size]), m_zero(0)
{
  for (long i = 0; i < m_size; i++)
    m_data[i] = Vector<T>(m_size);
}

template <typename T>
DenseMatrix<T>::DenseMatrix(const DenseMatrix<T>& copy)
: m_size(copy.GetSize()), m_data(new Vector<T>[m_size]), m_zero(copy.GetTolerance())
{
  for (long i = 0; i < m_size; i++)
    m_data[i] = Vector<T>(copy.m_data[i]);
}

template <typename T>
DenseMatrix<T>::DenseMatrix(const Matrix<T> * copy)
: m_size(copy->GetSize()), m_data(new Vector<T>[m_size]), m_zero(copy->GetTolerance())
{
  for (long i = 0; i < m_size; i++)
    m_data[i] = Vector<T>(m_size);

  for (long i = 0; i < m_size; i++)
    for (long j = 0; j < m_size; j++)
      (*this)(i, j, (*copy)(i, j));
}

template <typename T>
DenseMatrix<T>::~DenseMatrix()
{
  delete[] m_data;
} 

template <typename T>
void DenseMatrix<T>::SetTolerance(const T tol)
{
  if (tol < 0)
    throw TolErr("void DenseMatrix<T>::SetTolerance");
  m_zero = tol;
  for (long i = 0; i < m_size; i++)
    m_data[i].SetTolerance(m_zero);
}

template <typename T>
T DenseMatrix<T>::GetTolerance() const
{
  return m_zero;
}

template <typename T>
T DenseMatrix<T>::GetSize() const
{
  return m_size;
}

template <typename T>
Matrix<T>& DenseMatrix<T>::Transpose()
{
  for (long i = 0; i < m_size - 1; i++)
  {
    for (long j = i + 1; j < m_size; j++)
    {
      T tmp = (*this)(i, j);
      (*this)(i, j, (*this)(j, i));
      (*this)(j, i, tmp);
    }
  }
  return *this;
}

template <typename T>
bool DenseMatrix<T>::IsDiagDom() const
{
  for (long i = 0; i < m_size; i++)
  {
    T sum = 0;
    for (long j = 0; j < m_size; j++)
    {
      if (i != j)
        sum += abs((*this)(i, j));
    }
    if (sum > abs((*this)(i, i)))
      return false;
  }
  return true;
}

template <typename T>
T DenseMatrix<T>::operator()(const long row, const long col) const
{
  if (row < 0 || row >= m_size)
    throw SubscriptErr(row, "T DenseMatrix<T>::operator()");
  if (col < 0 || col >= m_size)
    throw SubscriptErr(col, "T DenseMatrix<T>::operator()");

  return m_data[row][col];
}

template <typename T>
Matrix<T>& DenseMatrix<T>::operator()(const long row, const long col, const T val)
{
  if (row < 0 || row >= m_size)
    throw SubscriptErr(row, "Matrix<T>& DenseMatrix<T>::operator()");
  if (col < 0 || col >= m_size)
    throw SubscriptErr(col, "Matrix<T>& DenseMatrix<T>::operator()");
  
  m_data[row][col] = val;
  return *this;
}

template <typename T>
unique_ptr<Matrix<T>> DenseMatrix<T>::operator*(const T rhs) const
{
  unique_ptr<Matrix<T>> m(new DenseMatrix<T>(*this));
  for (long i = 0; i < m_size; i++)
    for (long j = 0; j < m_size; j++)
      (*m)(i, j, (*m)(i, j) * rhs);
  return m;
}

template <typename T>
Vector<T> DenseMatrix<T>::operator*(const Vector<T>& rhs) const
{
  Vector<T> vect = Vector<T>(m_size);
  for (long i = 0; i < m_size; i++)
  {
    T sum = 0;
    for (long j = 0; j < m_size; j++)
    {
      sum += ((*this)(i, j) * rhs[j]);
    }
    vect[i] = sum;
  }
  return vect;
}

template <typename T>
unique_ptr<Matrix<T>> DenseMatrix<T>::operator*(const Matrix<T>& rhs) const
{
  unique_ptr<Matrix<T>> m(new DenseMatrix<T>(*this));
  for (long i = 0; i < m_size; i++)
  {
    for (long j = 0; j < m_size; j++)
    {
      T sum = 0;
      for (long k = 0; k < m_size; k++)
      {
        sum += ((*this)(i, k) * rhs(k, j));
      }
      (*m)(i, j, sum);
    }
  }
  return m;
}

template <typename T>
unique_ptr<Matrix<T>> DenseMatrix<T>::operator+(const T rhs) const
{
  //DenseMatrix<T> * m = new DenseMatrix<T>(*this);
  unique_ptr<Matrix<T>> m(new DenseMatrix<T>(*this));
  for (long i = 0; i < m_size; i++)
    for (long j = 0; j < m_size; j++)
      (*m)(i, j, (*m)(i, j) + rhs);
  return m;
}

template <typename T>
unique_ptr<Matrix<T>> DenseMatrix<T>::operator+(const Matrix<T>& rhs) const
{
  //DenseMatrix<T> * m = new DenseMatrix<T>(*this);
  unique_ptr<Matrix<T>> m(new DenseMatrix<T>(*this));
  for (long i = 0; i < m_size; i++)
    for (long j = 0; j < m_size; j++)
      (*m)(i, j, (*m)(i, j) + rhs(i, j));
  return m;
}

template <typename T>
Matrix<T>& DenseMatrix<T>::operator+=(const T rhs)
{
  *this = *(*this + rhs);
  return *this;
}

template <typename T>
Matrix<T>& DenseMatrix<T>::operator+=(const Matrix<T>& rhs)
{
  *this = *(*this + rhs);
  return *this;
}

template <typename T>
unique_ptr<Matrix<T>> DenseMatrix<T>::operator-(const T rhs) const
{
  unique_ptr<Matrix<T>> m(new DenseMatrix<T>(*this));
  *m += (rhs * -1);
  return m;
}

template <typename T>
unique_ptr<Matrix<T>> DenseMatrix<T>::operator-(const Matrix<T>& rhs) const
{
  unique_ptr<Matrix<T>> m(new DenseMatrix<T>(*this));
  *m += *(rhs * -1);
  return m;
}

template <typename T>
Matrix<T>& DenseMatrix<T>::operator-=(const T rhs)
{
  *this += (rhs * -1);
  return *this;
}

template <typename T>
Matrix<T>& DenseMatrix<T>::operator-=(const Matrix<T>& rhs)
{
  *this += *(rhs * -1);
  return *this;
}

template <typename T>
Matrix<T>& DenseMatrix<T>::operator=(const Matrix<T>& rhs)
{
  m_size = rhs.GetSize();
  m_zero = rhs.GetTolerance();
  for (long i = 0; i < m_size; i++)
    for (long j = 0; j < m_size; j++)
      (*this)(i, j, rhs(i, j));
  return *this;
}

template <typename T>
DenseMatrix<T>& DenseMatrix<T>::operator=(const DenseMatrix<T>& rhs)
{
  m_size = rhs.GetSize();
  m_zero = rhs.GetTolerance();
  for (long i = 0; i < m_size; i++)
    for (long j = 0; j < m_size; j++)
      (*this)(i, j, rhs(i, j));
  return *this;
}