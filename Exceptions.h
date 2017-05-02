#pragma once

#include <stdexcept>
#include <iostream>
#include <string>
using namespace std;

class Error
{
  public:
    virtual ~Error() {};
    virtual string what() const = 0;
};

class TolErr : public Error
{
  public:
    // Desc: default constructor for the TolErr class
    // Pre: None
    // Post: an instance of the TolErr class
    // will be created
    TolErr(const string w = "");

    // Desc: Default destructor
    // Pre: None
    // Post: The error will be deleted
    ~TolErr();

    // Desc: returns a string description of the error
    // Pre: None
    // Post: returns a string description of the error
    virtual string what() const;
  private:
    string where;
};

class DiagDomErr : public Error
{
  public:
    // Desc: default constructor for the DiagDomErr class
    // Pre: None
    // Post: an instance of the DiagDomErr class
    // will be created
    DiagDomErr(const string w = "");
    
    // Desc: Default destructor
    // Pre: None
    // Post: The error will be deleted
    ~DiagDomErr();

    // Desc: returns a string description of the error
    // Pre: None
    // Post: returns a string description of the error
    virtual string what() const;
  private:
    string where;
};

class SubscriptErr : public Error
{
  public:
    // Desc: default constructor for the SubscriptErr class
    // Pre: None
    // Post: an instance of the SubscriptErr class
    // will be created
    SubscriptErr(const long i = -1, const string w = "");
    
    // Desc: Default destructor
    // Pre: None
    // Post: The error will be deleted
    ~SubscriptErr();

    // Desc: returns a string description of the error
    // Pre: None
    // Post: returns a string description of the error
    virtual string what() const;
  private:
    long index;
    string where;
};

class DivByZeroErr : public Error
{
  public:
    // Desc: default constructor for DivByZeroErr class
    // Pre: None
    // Post: creates an instance of the DivByZeroErr class
    DivByZeroErr(const string w = "");
    
    // Desc: Default destructor
    // Pre: None
    // Post: The error will be deleted
    ~DivByZeroErr();

    // Desc: returns a string description of the error
    // Pre: None
    // Post: returns a string description of the error
    virtual string what() const;
  private:
    string where;
};

class SizeErr : public Error
{
  public:
    // Desc: default constructor for the SizeErr class
    // Pre: None
    // Post: an instance of the SizeErr class
    // will be created
    SizeErr(const long a = -1, const long b = -1, const string w = "");
    
    // Desc: Default destructor
    // Pre: None
    // Post: The error will be deleted
    ~SizeErr();

    // Desc: returns a string description of the error
    // Pre: None
    // Post: returns a string description of the error
    virtual string what() const;
  private:
    long m_a;
    long m_b;
    string where;
};
