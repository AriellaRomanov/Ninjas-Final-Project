/*********************************************
** Programmers: Andrea Smith, Brian Yadamec
** Class: Object-Oriented Numerical Modeling
** Project: Dirichlet Problem
** Date: May 7th, 2017
*********************************************/

#include <chrono>
#include "DenseMatrix.h"
#include "PDEMatrix.h"
#include "PDE.h"
using Type = double;

Type f_x0(Type a){return a*a*a;}
Type f_x1(Type a){return a*a*a;}
Type f_y0(Type a __attribute__((unused))){return 0.0;}
Type f_y1(Type a __attribute__((unused))){return 1.0;}
Type f_g(Type y){return y;}


int main(int argc, char * argv[])
{  
  if (argc < 3)
  {
    cout << "Program usage: " << argv[0] << " <N> <zero tolerance> <output file (optional)>" << endl;
    return 1;
  }

  try
  {
    long N = atoi(argv[1]);
    if (N < 1)
      throw SizeErr(N, 1, "Expected min N of 1");
    Type tol = atof(argv[2]);
    if (tol < 0)
      throw TolErr("program input");
  
    PDESolution<Type> X(f_x0, f_x1, f_y0, f_y1, f_g);
    auto timeStart = chrono::high_resolution_clock::now();
    Vector<Type> jvect(X.Jacobi(N, tol));
    auto jDiff = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - timeStart).count();

    timeStart = chrono::high_resolution_clock::now();
    Vector<Type> gvect(X.Gaussian(N, tol));
    auto gDiff = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - timeStart).count();

    if (argc >= 4)
    {
      Vector<Type> x(jvect.GetSize());
      Vector<Type> y(jvect.GetSize());
      for (long i = 0; i < N - 1; i++)
        for (long j = 0; j < N - 1; j++)
        {
          x[i * (N - 1) + j] = (j + 1) / static_cast<Type>(N);
          y[i * (N - 1) + j] = (i + 1) / static_cast<Type>(N);
        }

      //file
      ofstream file;
      file.open(argv[3]);
      if (file.is_open())
      {
        file << x << endl;
        file << y << endl;
        file << jvect << endl << endl;
        file << gvect << endl;
        file << "Jacobi (" << jDiff << " milliseconds)" << endl;
        file << "Gaussian (" << gDiff << " milliseconds)" << endl;
        file.close();
      }
      else
        cout << "Could not open file " << argv[3] << endl;
    }
    else
    {
      cout << jvect << endl << endl;
      cout << gvect << endl;
    }
    cout << "Jacobi (" << jDiff << " milliseconds)" << endl;
    cout << "Gaussian (" << gDiff << " milliseconds)" << endl;
  
  }
  catch (Error& err)
  {
    cout << err.what() << endl;
  }

  return 0;
}