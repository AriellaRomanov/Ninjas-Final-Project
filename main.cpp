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
    long timeStart = time(NULL);
    Vector<Type> jvect(X.Jacobi(N, tol));
    long jDiff = time(NULL) - timeStart;

    timeStart = time(NULL);
    Vector<Type> gvect(X.Gaussian(N, tol));
    long gDiff = time(NULL) - timeStart;
    
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
        //file << "Jacobi (" << jDiff << " seconds):" << endl;
        
        file << x << endl;
        file << y << endl;
        file << jvect << endl << endl;

        //for (long i = 0; i < (N - 1) * (N - 1); i++)
          //file << x[i] << " " << y[i] << " " << jvect[i] << endl;

        // file << "ListPointPlot3D[{";
        // for (long i = 0; i < (N - 1) * (N - 1); i++)
        // {
        //   file << "{" << x[i] << "," << y[i] << "," << jvect[i] << "}";
        //   if (i < (N - 1) * (N - 1) - 1)
        //     file << ",";
        // }
        // file << "}]" << endl;

        //file << "Gaussian (" << gDiff << " seconds):" << endl;

        file << gvect << endl;

        //for (long i = 0; i < (N - 1) * (N - 1); i++)
          //file << x[i] << " " << y[i] << " " << gvect[i] << endl;

        // file << "ListPointPlot3D[{";
        // for (long i = 0; i < (N - 1) * (N - 1); i++)
        // {
        //   file << "{" << x[i] << "," << y[i] << "," << gvect[i] << "}";
        //   if (i < (N - 1) * (N - 1) - 1)
        //     file << ",";
        // }
        // file << "}]" << endl;
        
        file.close();
      }
      else
        cout << "Could not open file " << argv[3] << endl;
    }
    else
    {
      cout << "Jacobi (" << jDiff << " seconds):" << endl;
      cout << jvect << endl;
      cout << "Gaussian (" << gDiff << " seconds):" << endl;
      cout << gvect << endl;
    }
  
  }
  catch (Error& err)
  {
    cout << err.what() << endl;
  }

  return 0;
}