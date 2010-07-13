// An interface for different local optimisation algorithms.
// The main purpose of the module is to hide some global variables
// we have to define in order to use external optimisation functions.

#ifndef __LOCAL_HH
#define __LOCAL_HH

#include "common.hh"
#include "coating.hh"
#include "target.hh"
#include <vector>


void ConjugateGradientOptimisation(Coating& coating,
				   AbstractTarget& target);

void ConjugateGradientOptimisation(vector<Coating>&,
				   AbstractTarget& target);

void SimplexOptimisation(Coating& coating, AbstractTarget& target, Real step);
/* step is the value we add to the thicknesses of the original design
   in order to form the initial simplex */

void SimplexOptimisation(vector<Coating>& coatings,
			 AbstractTarget& target, Real step);

#endif
