#ifndef __RANDOM_HH
#define __RANDOM_HH

#include "common.hh"


unsigned int InitialiseRandomNumberGenerator(unsigned int seed = 0);
// initialises the random number generator; note that
// if you call any of the following functions without
// calling first this one, it will be called automatically;
// if the seed is equal to zero, it uses the current time
// to initialise the generator; returns the seed used

Real UniformRandom(Real min = Real(0), Real max = Real(1));
// returns a uniformly distributed random number between
// min and max

int RandomInteger(int min, int max);
// returns a uniformly distributed random number between
// min and max

int UnbiasedCoin();
//  returns either 1 or -1 with equal probabilities

Real GaussianRandom();
// returns a normally distributed deviate with zero mean and
// unit variance

void Permutation(int N, int* array);
// makes a random permutation of integers from 0 to (N-1) and
// writes it into the 'array'

#endif
