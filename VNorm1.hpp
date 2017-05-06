VNorm1::VNorm1() {}

template <typename T>
T VNorm1::operator()(const Vector<T>& vect)
{
  T result = 0;
  for (long i = 0; i < vect.GetSize(); i++)
    result += abs(vect[i]);
  return result;
}