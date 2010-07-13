#ifndef __DISPERSION_HH
#define __DISPERSION_HH

#include "common.hh"
#include <vector>

class AbstractDispersion
{
public:
    AbstractDispersion() {}
    virtual ~AbstractDispersion() {}
    virtual Real Phase(Real omega) const = 0;
    virtual Real GD(Real omega) const = 0;
    virtual Real GDD(Real omega) const = 0;
};

//-------------------------------------------------------------------------

class TaylorDispersion : public AbstractDispersion
{
private:
    Real omega0;
    vector<Real> beta; // the first element in this array is
                       // the GDD at omega0, the second
                       // element is the TOD and so on
public:
    TaylorDispersion();
    TaylorDispersion(Real omega0, vector<Real> coefficients);
    ~TaylorDispersion();
    void Set(Real omega0, vector<Real> coefficients);
    Real Omega0() { return omega0; }
    virtual Real Phase(Real omega) const;
    virtual Real GD(Real omega) const;
    virtual Real GDD(Real omega) const;
};

//-------------------------------------------------------------------------

Real AdaptiveDispersionFactor(int N, const Real* target_GDD,
			      const Real* GDD_reserve,
			      const Real* GDD_tolerance,
			      const Real* GDD);
/* returns the factor, by which the target dispersion has to be multiplied
   in order to match the given dispersion */

Real GDShift(int N, const Real* target_GD, const Real* GD_reserve,
	      const Real* GD_tolerance, const Real* GD);
/* returns the value, which has to be added to 'GD' in order to make it
   as close as possible to the target group delay 'target_GD' */


#endif
