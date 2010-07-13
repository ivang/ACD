#ifndef __COMMON_H
#define __COMMON_H

#include <complex>
using namespace std;

#ifdef DEBUG
#include <cstdlib>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#endif

#ifdef MS_WINDOWS
#define DEFAULT_PATH_TO_MATERIALS "C:\\ACD\\MATERIALS\\"
const char directory_separator = '\\';
#else
#define DEFAULT_PATH_TO_MATERIALS "/usr/local/share/acd/materials/"
const char directory_separator = '/';
#endif

enum Tpolarisation {TE, TM};

#ifdef SINGLE_PRECISION
typedef float Real;
typedef complex<float> Complex;
#else
typedef double Real;
typedef complex<double> Complex;
#endif

const Real SPEED_OF_LIGHT=299.792458;

//-------------------------------------------------------------------------

template <class T> inline T square(T x) { return x*x; }

template <class T>
inline T IntPower(T x, int n) // returns x^n
{
	T result;
	if (x == T(0)) return 0;
	if (n == 0) return static_cast<T>(1);
	else if (n < 0)
	{
		x = static_cast<T>(1) / x;
		n = - n;
	}
	result = static_cast<T>(1);
	while (n > 0)
	{
		if (n % 2)
		{
			result *= x;
			n--;
		}
		else
		{
			x *= x;
			n /= 2;
		}
	}
	return result;
}

//-------------------------------------------------------------------------

#endif
