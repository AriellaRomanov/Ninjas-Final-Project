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

template <typename T>
Vector<T> DenseMatrix<T>::Gaussian(long m_n, Vector<T>& b_vector)
{
  if(m_size != (m_n-1)*(m_n-1))
    throw SizeErr(m_size,(m_n-1)*(m_n-1),"Gaussian Matrix/Passed Size");
  if(m_size != b_vector.GetSize())
    throw SizeErr(m_size,b_vector.GetSize(),"Gaussian Matrix/Vector");
  for(long row=0;row<(m_n-1)*(m_n-1);row++)
  {
    for(long col=0;col<(m_n-1)*(m_n-1);col++)
    {
      m_data[row][col] =  0;
      //calculate diagonal
      if (row == col)
        m_data[row][col] = 1;

      //calculate lower off diagonal
      if (row == col + 1)
      {
      if (row % (m_n - 1) == 0)
        m_data[row][col] =  0;
      else
        m_data[row][col] =  -0.25;
      }

      //calculate upper off diagonal
      if (col == row + 1)
      {
        if (col % (m_n - 1) == 0)
          m_data[row][col] =  0;
        else
          m_data[row][col] =  -0.25;
      }

      //further diagonals
      if (row == col + m_n - 1 || col == row + m_n - 1)
        m_data[row][col] =  -0.25;      
    }
  }
  
  
  
  for(long i=0;i<(m_n-1)*(m_n-1);i++)
  {
    for(long j=i+1;j<(m_n-1)*(m_n-1);j++)
    {
      if(m_data[i][j] != 0)
      {
        T temp = m_data[i][j]/m_data[i][i];
        for(long k=i;k<(m_n-1)*(m_n-1);k++)
          m_data[k][j] -= temp*m_data[k][i];
        b_vector[j] -= temp*b_vector[i];
      }
    }
  }
  //Bottom half is 0's
  for(long i=((m_n-1)*(m_n-1))-1;i>=0;i--)
  {
    for(long k=0;k<i;k++)
    {
      b_vector[k] -= (m_data[i][k]/m_data[i][i])*b_vector[i];
      m_data[i][k] = 0;
    }
  }
  //Top half is 0's
  for(long i=0;i<(m_n-1)*(m_n-1);i++)
  {
    if(m_data[i][i] == 0)
      throw  DivByZeroErr("Gaussian");
    b_vector[i] = b_vector[i]/m_data[i][i];
  }
  //Diagonal is 1's
  return b_vector;
}