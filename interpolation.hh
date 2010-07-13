#ifndef __INTERPOLATION_HH
#define __INTERPOLATION_HH

#include "common.hh"
#include <vector>

class AbstractInterpolator
{
protected:
    int n;
    Real *X_copy; // a copy of 'X'
    Complex *Y_copy; // a copy of 'Y'
    void Initialise(const vector<Real>& X, const vector<Complex>& Y);
    // 'Initialise' allocates the memory for the arrays
    // and copies the data there
public:
    AbstractInterpolator();
    AbstractInterpolator(const vector<Real>& X, const vector<Complex>& Y);
    virtual ~AbstractInterpolator();
    virtual void Reset(const vector<Real>& X, const vector<Complex>& Y);
    virtual Complex operator() (Real x) = 0;
};


class LinearInterpolator : public AbstractInterpolator
{
public:
    LinearInterpolator();
    LinearInterpolator(const vector<Real>& X, const vector<Complex>& Y);
    ~LinearInterpolator();
    void Reset(const vector<Real>& X, const vector<Complex>& Y);
    Complex operator() (Real x);
};


class SplineInterpolator : public AbstractInterpolator
{
private:
    Complex* D2; // the second derivative y''(x)
    void Initialise(const vector<Real>& X, const vector<Complex>& Y);
public:
    SplineInterpolator();
    SplineInterpolator(const vector<Real>& X, const vector<Complex>& Y);
    ~SplineInterpolator();
    void Reset(const vector<Real>& X, const vector<Complex>& Y);
    Complex operator() (Real x);
};

#endif
