#include "dispersion.hh"
// #include "reader.hh"
// #include <cmath>


TaylorDispersion::TaylorDispersion(Real omega0,
				   vector<Real> coefficients) :
    AbstractDispersion(),
    omega0(omega0), beta(coefficients)
{
}

TaylorDispersion::TaylorDispersion() : AbstractDispersion()
{
    omega0 = 0.0;
}

void TaylorDispersion::Set(Real _omega0, vector<Real> coefficients)
{
    omega0 = _omega0;
    beta = coefficients;
}

TaylorDispersion::~TaylorDispersion() {}

Real TaylorDispersion::Phase(Real omega) const
{
    int i, n;
    Real x = omega - omega0;
    Real phi = 0;
    n = beta.size();
    for (i=n-1; i>=0; i--) phi = (phi + beta[i]) * x / Real(i+2);
    phi *= -x;
    return phi;
}

Real TaylorDispersion::GD(Real omega) const
{
    int i, n;
    Real x = omega - omega0;
    Real GD = 0;
    n = beta.size();
    for (i=n-1; i>=0; i--) GD = (GD + beta[i]) * x / Real(i+1);
    return GD;
}

Real TaylorDispersion::GDD(Real omega) const
{
    int i, n;
    Real x = omega - omega0;
    Real GDD = 0;
    n = beta.size();
    for (i=n-1; i>0; i--) GDD = (GDD + beta[i]) * x / Real(i);
    GDD += beta[0];
    return GDD;
}

//-------------------------------------------------------------------------

Real AdaptiveDispersionFactor(int N, const Real* target_GDD,
			      const Real* GDD_reserve,
			      const Real* GDD_tolerance, const Real* GDD)
{
    int i;
    Real x, s1=0, s2=0;
    for (i=0; i<N; i++)
    {
	if (abs(GDD[i] - target_GDD[i]) > GDD_reserve[i])
	{
	    x = square(GDD_tolerance[i]);
	    s1 += target_GDD[i] * GDD[i] / x;
	    s2 += square(target_GDD[i]) / x;
	}
    }
    if (s2==Real(0)) return 1;
    return s1/s2;
}


Real GDShift(int N, const Real* target_GD, const Real* GD_reserve,
	      const Real* GD_tolerance, const Real* GD)
{
    int i;
    Real x, s1=0, s2=0;
    for (i=0; i<N; i++)
    {
	if (abs(GD[i] - target_GD[i]) > GD_reserve[i])
	{
	    x = square(GD_tolerance[i]);
	    s1 += (target_GD[i] - GD[i]) / x;
	    s2 += 1 / x;
	}
    }
    if (s2==Real(0)) return 0;
    return s1/s2;
}

//=========================================================================
