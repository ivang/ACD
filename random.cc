#include "random.hh"
#include <cstdlib>
#include <cmath>
#include <time.h>

static bool generator_initialised = false;

unsigned int InitialiseRandomNumberGenerator(unsigned int seed)
{
    time_t t;
    if (! seed) seed = (unsigned int) time(&t);
#ifdef MS_WINDOWS
    srand(seed);
#else
    srandom(seed);
#endif
    generator_initialised = true;
    return seed;
}

Real UniformRandom(Real min, Real max)
{
    if (! generator_initialised)
    {
	InitialiseRandomNumberGenerator();
    }
#ifdef MS_WINDOWS
    return min + (Real)rand() / (RAND_MAX + 1.0) * (max - min);
#else
    return min + (Real)random() / (RAND_MAX + 1.0) * (max - min);
#endif
}

int RandomInteger(int min, int max)
{
    if (! generator_initialised)
    {
	InitialiseRandomNumberGenerator();
    }
#ifdef MS_WINDOWS
    return min + rand() % (max - min + 1);
#else
    return min + random() % (max - min + 1);
#endif
}

Real GaussianRandom()
{	/* returns a normally distributed deviate with zero mean and unit
	   variance */
    static int iset = 0;
    static Real gset;
    Real fac,rsq,v1,v2;
    if (! generator_initialised)
    {
	InitialiseRandomNumberGenerator();
    }
    if (iset == 0)
    {
	do
	{
#ifdef MS_WINDOWS
	    v1=2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
	    v2=2.0 * rand() / (RAND_MAX + 1.0) - 1.0;
#else
	    v1=2.0 * random() / (RAND_MAX + 1.0) - 1.0;
	    v2=2.0 * random() / (RAND_MAX + 1.0) - 1.0;
#endif
	    rsq=v1*v1+v2*v2;
	} while (rsq >= Real(1) || rsq == Real(0));
	fac = sqrt(-2.0*log(rsq)/rsq);
	gset = v1*fac;
	iset = 1;
	return v2*fac;
    }
    else
    {
	iset=0;
	return gset;
    }
}

//-------------------------------------

void Permutation(int N, int* array)
{  /* makes a random permutation of integers from 0 to (N-1) and
      writes it into the 'array' */
    int i, j, swap;
    if (! generator_initialised)
    {
	InitialiseRandomNumberGenerator();
    }
    for (i = 0; i < N; i++) array[i] = i;
    for (i = N-1; i >= 0; i--)
    {
#ifdef MS_WINDOWS
	j = rand() % (i+1);
#else
	j = random() % (i+1);
#endif
	swap = array[j];
	array[j] = array[i];
	array[i] = swap;
    }
}


int UnbiasedCoin()
{
#ifdef MS_WINDOWS
    if (rand() > RAND_MAX/2) return 1;
#else
    if (random() > RAND_MAX/2) return 1;
#endif
    return -1;
}
