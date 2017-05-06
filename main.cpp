#include "DenseMatrix.h"
#include "PDEMatrix.h"
#include "PDE.h"
using Type = double;

/* template <typename T, T T_funct(T)>
T evaluate(T operand) 
{
  return T_funct(operand);
}
 */
Type f_x0(Type a){return a*a*a;}
Type f_x1(Type a){return a*a*a;}
Type f_y0(Type a __attribute__((unused))){return 0.0;}
Type f_y1(Type a __attribute__((unused))){return 1.0;}
Type f_g(Type y){return y;}


int main()
{
  try
  {
    long size = 10;
    long iteration_max = 1000000;
    
    PDESolution<Type> X(f_x0, f_x1, f_y0, f_y1, f_g);
    int timeStart = time(NULL);
    cout << X.Jacobi(iteration_max, size) << endl;
    cout<<time(NULL)-timeStart<<endl;
    
    timeStart = time(NULL);
    cout << X.Gaussian(size) << endl;
    cout<<time(NULL)-timeStart<<endl;
    
    /* const long SIZE = 4;
    Vector<Type> b_vector((SIZE-1)*(SIZE-1));
    for(int i=0;i<(SIZE-1)*(SIZE-1);i++)
      b_vector[i] = 0;
    for(int i=0;i<SIZE-1;i++)
      b_vector[i] += (evaluate<Type,x0>(static_cast<double>((i+1))/(SIZE)))/4;
    for(int i=0;i<SIZE-1;i++)
      b_vector[i*(SIZE-1)] += (evaluate<Type,y0>(static_cast<double>((i+1))/(SIZE)))/4;
    for(int i=0;i<SIZE-1;i++)
      b_vector[(((SIZE-1)*(SIZE-1))-(SIZE-1))+i] += (evaluate<Type,x1>(static_cast<double>((i+1))/(SIZE)))/4;
    for(int i=0;i<SIZE-1;i++)
      b_vector[((i+1)*(SIZE-1))-1] += (evaluate<Type,y1>(static_cast<double>((i+1))/(SIZE)))/4;
    
    cout<<b_vector<<endl; */
  
  }
  catch (Error& err)
  {
    cout << err.what() << endl;
  }

  return 0;
}