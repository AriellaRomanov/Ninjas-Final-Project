#include "DenseMatrix.h"
using Type = double;

int main()
{
  const long SIZE = 5;

  try
  {
    DenseMatrix<Type> matrix(SIZE);
    for (long i = 0; i < SIZE; i++)
      for (long j = 0; j < SIZE; j++)
        matrix(i, j, i * j);

    cout << matrix(SIZE, SIZE) << endl;
  }
  catch (Error& err)
  {
    cout << err.what() << endl;
  }

  return 0;
}