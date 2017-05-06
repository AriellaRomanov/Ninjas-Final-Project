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
Vector<T> PDESolution<T>::Gaussian(long m_n)
{
  long m_size = m_n-1;
  Vector<T> b_vector((m_n-1)*(m_n-1));
  for(int i=0;i<(m_n-1)*(m_n-1);i++)
    b_vector[i] = 0;
  for(int i=0;i<m_n-1;i++)
    b_vector[i] += m_x0(static_cast<T>(i+1)/(m_n))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[i*(m_n-1)] += m_y0(static_cast<T>(i+1)/(m_n))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[(((m_n-1)*(m_n-1))-(m_n-1))+i] += m_x1(static_cast<T>(i+1)/(m_n))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[((i+1)*(m_n-1))-1] += m_y1(static_cast<T>(i+1)/(m_n))/4;
  for(int i=0;i<m_size;i++)
  {
    for(int k=0;k<m_size;k++)
      b_vector[(i*m_size)+k] -= m_g(static_cast<T>(i+1)/m_n)*(1/static_cast<T>(m_n*m_n*4));
  }
  
  DenseMatrix<T> A(m_size*m_size);
  return A.Gaussian(m_n,b_vector);
}

template <typename T>
Vector<T> PDESolution<T>::Jacobi(long iter_max, long m_n)
{
  long m_size = m_n-1;
  Vector<T> b_vector((m_n-1)*(m_n-1));
  for(int i=0;i<(m_n-1)*(m_n-1);i++)
    b_vector[i] = 0;
  for(int i=0;i<m_n-1;i++)
    b_vector[i] += m_x0(static_cast<T>(i+1)/(m_n))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[i*(m_n-1)] += m_y0(static_cast<T>(i+1)/(m_n))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[(((m_n-1)*(m_n-1))-(m_n-1))+i] += m_x1(static_cast<T>(i+1)/(m_n))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[((i+1)*(m_n-1))-1] += m_y1(static_cast<T>(i+1)/(m_n))/4;
  for(int i=0;i<m_size;i++)
  {
    for(int k=0;k<m_size;k++)
      b_vector[(i*m_size)+k] -= m_g(static_cast<T>(i+1)/m_n)*(1/static_cast<T>(m_n*m_n*4));
  }
  
  Vector<T> ErrT(m_size*m_size);
  for(int i=0;i<m_size*m_size;i++)
    ErrT[i] = .0000001;
	Vector<T> x_prev(m_size*m_size);
  PDEMatrix<T> a_matrix(m_n);
  
  for(int i=0;i<iter_max;i++)
	{
	  Vector<T> temp(m_size*m_size);
	  temp = b_vector - a_matrix.Jacob_Mult(x_prev);
    x_prev = temp;
    
	  if(((a_matrix*temp) - b_vector) < ErrT)
	  {
	    cout<<"Solution after "<<i<<" iterations."<<endl;
      i=iter_max;
	  }
	}
  
  return x_prev;
}