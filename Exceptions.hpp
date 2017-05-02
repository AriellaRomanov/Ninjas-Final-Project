TolErr::TolErr(const string w)
: where(w)
{}

TolErr::~TolErr()
{};

string TolErr::what() const
{
  return "Cannot set a negative tolerance at " + where;
}
    
DiagDomErr::DiagDomErr(const string w)
: where(w)
{}

DiagDomErr::~DiagDomErr()
{};

string DiagDomErr::what() const
{
  return "Matrix not diagonal dominant at " + where;
}

SubscriptErr::SubscriptErr(const long i, const string w)
: index(i), where(w)
{}

SubscriptErr::~SubscriptErr()
{};

string SubscriptErr::what() const
{
  string err = "Bad subscript: " +  to_string(index) + " at " + where;
  return err;
}

DivByZeroErr::DivByZeroErr(const string w)
: where(w)
{}

DivByZeroErr::~DivByZeroErr()
{};

string DivByZeroErr::what() const
{
  return "Illegal division by zero at " + where;
}

SizeErr::SizeErr(const long a, const long b, const string w)
: m_a(a), m_b(b), where(w)
{}

SizeErr::~SizeErr()
{};

string SizeErr::what() const
{
  string err = "Non-matching sizes: " +  to_string(m_a) + " and " + to_string(m_b) + " at " + where;
  return err;
}