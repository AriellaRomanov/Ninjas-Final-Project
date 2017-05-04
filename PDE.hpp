
template <typename T>
PDESolution::PDESolution(long n, T x0(T), T x1(T), T y0(T), T y1(T), T g(T)
{
  m_n = n;
  m_size = n-1;
  m_x0 = x0;
  m_x1 = x1;
  m_y0 = y0;
  m_y1 = y1;
  m_g = g; 
}

template <typename T>
Vector<T> PDESolution::Guassian()
{
  Vector<Type> b_vector((m_n-1)*(m_n-1));
  for(int i=0;i<(m_n-1)*(m_n-1);i++)
    b_vector[i] = 0;
  for(int i=0;i<m_n-1;i++)
    b_vector[i] += (evaluate<Type,x0>(static_cast<double>((i+1))/(m_n)))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[i*(m_n-1)] += (evaluate<Type,y0>(static_cast<double>((i+1))/(m_n)))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[(((m_n-1)*(m_n-1))-(m_n-1))+i] += (evaluate<Type,x1>(static_cast<double>((i+1))/(m_n)))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[((i+1)*(m_n-1))-1] += (evaluate<Type,y1>(static_cast<double>((i+1))/(m_n)))/4;
  //Not finished
  
}

template <typename T>
Vector<T> PDESolution::Jacobi(long iter_max)
{
  Vector<T> b_vector((m_n-1)*(m_n-1));
  for(int i=0;i<(m_n-1)*(m_n-1);i++)
    b_vector[i] = 0;
  for(int i=0;i<m_n-1;i++)
    b_vector[i] += (evaluate<T,x0>(static_cast<T>((i+1))/(m_n)))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[i*(m_n-1)] += (evaluate<T,y0>(static_cast<T>((i+1))/(m_n)))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[(((m_n-1)*(m_n-1))-(m_n-1))+i] += (evaluate<T,x1>(static_cast<T>((i+1))/(m_n)))/4;
  for(int i=0;i<m_n-1;i++)
    b_vector[((i+1)*(m_n-1))-1] += (evaluate<T,y1>(static_cast<T>((i+1))/(m_n)))/4;
  
  
  //Vector<T> ErrT(A.GetSize(), ErrTolerance);
  //Error Tolerance needs to be set up
	Vector<T> x_prev(m_size);
  PDEMatrix a_matrix(m_n);
  for(int i=0;i<iter_max;i++)
	{
	  my_vector<float> temp(m_size);
	  temp = b_vector - a_matrix.Jacob_Mult(x_prev);
	  /* if(((A*temp) - b) < ErrT)
	  {
	    cout<<"Solution after "<<i<<" iterations."<<endl;
		i=iter_max;
	  } */
	  x_prev = temp;
	}
}

//Not sure about this templating
template <typename T, T T_funct(T)>
T PDESolution::evaluate(T operand) 
{
  return T_funct(operand);
}