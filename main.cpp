#include "DenseMatrix.h"
#include "PDEMatrix.h"
using Type = double;

template <typename T, T T_funct(T)>
T evaluate(T operand)
{
  return T_funct(operand);
}

Type x0(Type a){return a*a*a;}
Type x1(Type a){return a*a*a;}
Type y0(){return 0.0;}
Type y1(){return 1.0;}
Type g(Type y){return y;}


int main()
{
  try
  {
    
    const long SIZE = 4;
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
    
    //cout<<b_vector<<endl;
  
  }
  catch (Error& err)
  {
    cout << err.what() << endl;
  }

  return 0;
}