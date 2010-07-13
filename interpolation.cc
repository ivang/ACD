#include "interpolation.hh"

AbstractInterpolator::AbstractInterpolator()
{
    n = 0;
    X_copy = NULL;
    Y_copy = NULL;
};

AbstractInterpolator::
AbstractInterpolator(const vector<Real>& X, const vector<Complex>& Y)
{
    Initialise(X, Y);
}

AbstractInterpolator::~AbstractInterpolator()
{
    if (n>0)
    {
	delete[] X_copy;
	delete[] Y_copy;
    }
}

void AbstractInterpolator::
Initialise(const vector<Real>& X, const vector<Complex>& Y)
{
    int i;
    n = X.size();
    if (n != (int)Y.size())
	throw "ERROR in AbstractInterpolator::Initialise: \
X and Y have different sizes";
    for (i=1; i<n; i++)
    {
	if (X[i-1] >= X[i]) throw "ERROR in AbstractInterpolator::Initialise: \
'X' must be sorted in the ascending order";
    }
    X_copy = new Real[n];
    Y_copy = new Complex[n];
    for (i=0; i<n; i++)
    {
	X_copy[i] = X[i];
	Y_copy[i] = Y[i];
    }
}

void AbstractInterpolator::
Reset(const vector<Real>& X, const vector<Complex>& Y)
{
    if (n != 0)
    {
	delete[] X_copy;
	delete[] Y_copy;
    }
    Initialise(X, Y);
}

//-------------------------------------------------------------------------

LinearInterpolator::LinearInterpolator() : AbstractInterpolator() {}

LinearInterpolator::
LinearInterpolator(const vector<Real>& X, const vector<Complex>& Y) :
    AbstractInterpolator(X, Y) {}

LinearInterpolator::~LinearInterpolator() {}

void LinearInterpolator::
Reset(const vector<Real>& X, const vector<Complex>& Y)
{
    AbstractInterpolator::Reset(X, Y);
}

Complex LinearInterpolator::operator() (Real x)
{
    Real h;
    int k, k1, k2;
    if (n==1) return Y_copy[0];
    if (x<X_copy[0]) return Y_copy[0];
    if (x>X_copy[n-1]) return Y_copy[n-1];
    k1 = 0;
    k2 = n-1;
    // we will find the right place in the table by means of bisection
    while (k2-k1 > 1)
    {
	k = (k2+k1) / 2;
	if (X_copy[k] > x) k2 = k;
	else k1 = k;
    }
    // k1 and k2 now bracket the input value of x.
    h = X_copy[k2] - X_copy[k1];
    if (h == 0.) throw "ERROR in LinearInterpolator: \
the X's must be distinct";
    return Y_copy[k1] + (Y_copy[k2]-Y_copy[k1])* (x-X_copy[k1])/h;
}


//-------------------------------------------------------------------------

SplineInterpolator::SplineInterpolator() : AbstractInterpolator()
{
    D2 = NULL;
}

SplineInterpolator::SplineInterpolator
(const vector<Real>& X, const vector<Complex>& Y) :
    AbstractInterpolator(X, Y)
{
    Initialise(X, Y);
}

SplineInterpolator::~SplineInterpolator()
{
    if (n>0) delete[] D2;
}

void SplineInterpolator::
Initialise(const vector<Real>& X, const vector<Complex>& Y)
{ // we assume the "natural spline" here
    int i, k;
    Real sig;
    Complex p;
    Complex *u;

    D2 = new Complex[n];
    u = new Complex[n];
    D2[0] = 0.0;
    u[0] = 0.0;
    for (i=1; i<n-1; i++)
    {
	// This is the decomposition loop of the tridiagonal algorithm.
	// D2 and u are used for temporary storage of the decomposed factors. 
	sig = (X[i]-X[i-1]) / (X[i+1]-X[i-1]);
	p = sig * D2[i-1] + Real(2);
	D2[i] = (sig - Real(1)) / p;
	u[i] = (Real(6)*((Y[i+1]-Y[i]) / (X[i+1]-X[i]) - 
		     (Y[i]-Y[i-1]) / (X[i]-X[i-1])) / (X[i+1]-X[i-1]) -
		sig*u[i-1]) / p;
    }
    D2[n-1] = 0.0;
    for (k=n-2; k>=0; k--)
    {
	// This is the backsubstitution loop of the tridiagonal algorithm.
	D2[k] = D2[k]*D2[k+1] + u[k];
    }
    delete[] u;
}


void SplineInterpolator::
Reset(const vector<Real>& X, const vector<Complex>& Y)
{
    if (n>0) delete[] D2;
    AbstractInterpolator::Reset(X, Y);
    Initialise(X, Y);
}


Complex SplineInterpolator::operator() (Real x)
{
    Real a, b, h;
    int k, k1, k2;
    if (x<X_copy[0]) return Y_copy[0];
    if (x>X_copy[n-1]) return Y_copy[n-1];
    k1 = 0;
    k2 = n-1;
    // We will find the right place in the table by means of bisection. This
    // is optimal if sequential calls to this routine are at random values of
    // x. If sequential calls are in order, and closely spaced, one would do
    // better to store previous values of k1 and k2 and test if they remain
    // appropriate on the next call.
    while (k2-k1 > 1)
    {
	k = (k2+k1) / 2;
	if (X_copy[k] > x) k2 = k;
	else k1 = k;
    }
    // k1 and k2 now bracket the input value of x.
    h = X_copy[k2] - X_copy[k1];
    if (h == 0.) throw "ERROR in SplineInterpolator: \
the X's must be distinct";
    a = (X_copy[k2] - x) / h;
    b = (x - X_copy[k1]) / h; // Cubic spline polynomial is now evaluated.
    return a*Y_copy[k1]+b*Y_copy[k2] +
	(a*(a*a-1)*D2[k1]+b*(b*b-1)*D2[k2])*(h*h)/Real(6);
}
