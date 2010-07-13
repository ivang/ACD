#include "pulse.hh"
#include <cmath>

//------------------ general purpose functions ----------------------------

// Simpson integrator
template<class T>
T IntegrateArray(int n, T* A, Real dx)
{
  int i;
  T s1, s2;
  s1 = 0;
  s2 = 0;
  if (n % 2 == 1)
  {
    for (i=1; i<n-1; i += 2) s1 += A[i];
    for (i=2; i<n-1; i += 2) s2 += A[i];
  }
  else
  { // assume the function almost constant at the farther wing
    for (i=1; i<n; i += 2) s1 += A[i];
    for (i=2; i<n; i += 2) s2 += A[i];
  }
  return (A[0]+static_cast<T>(4)*s1+static_cast<T>(2)*s2+A[n-1])*dx/
      static_cast<T>(3);
}

void SmoothenPhase(Real* phase, int N)
{
    /* "Phase" must be a 1D array of real numbers. The procedure adds
       2*pi*<some integer number> to each element of the array passed
       in order to make it as smooth as possible; "N" is the number
       of elements in "phase". */
    Real x;
    int i, i0, k;
    if (N==1) return;
    i0 = N/2; // we'll try to make the phase at this point as small as possible
    if (i0+1 < N)
    {
	// correct the first point to the right
	x = phase[i0+1] - phase[i0];
	k = int(round(x/(2*M_PI)));
	phase[i0+1] -= 2*M_PI*k;
    }
    // correct the rest of points to the right from i0
    for (i=i0+2;i<N;i++)
    {
	x = phase[i] - 2.0*phase[i-1] + phase[i-2];
	k = int(round(x/(2*M_PI)));
	phase[i] -= 2*M_PI*k;
    }
    if (i0-1 >= 0)
    {
	// correct the first point to the left
	x = phase[i0-1] - phase[i0];
	k = int(round(x/(2*M_PI)));
	phase[i0-1] -= 2*M_PI*k;
    }
    // correct the rest of points to the left from i0
    for (i=i0-2; i>=0; i--)
    {
	x = phase[i] - 2.0*phase[i+1] + phase[i+2];
	k = int(round(x/(2*M_PI)));
	phase[i] -= 2*M_PI*k;
    }
}

//-------------------------------------------------------------------------

AbstractPulse::AbstractPulse(int Nt, Real t_step, bool in_time_domain) :
    Nt(Nt), t_step(t_step), is_in_time_domain(in_time_domain) {}

AbstractPulse::~AbstractPulse() {}

int AbstractPulse::NumberOfPoints() const { return Nt; }

Real AbstractPulse::TimeStep() const { return t_step; }

Real AbstractPulse::OmegaStep() const
{
    return 2.0*M_PI / (t_step*Nt);
}

bool AbstractPulse::IsInTimeDomain() const
{
    return is_in_time_domain;
}

bool AbstractPulse::IsInFrequencyDomain() const
{
    return (! is_in_time_domain);
}

//-------------------------------------------------------------------------


SimplePulse::SimplePulse(int Nt, Real t_step,
			 bool in_time_domain, bool fft_estimate) :
    AbstractPulse(Nt, t_step, in_time_domain), fft_estimate(fft_estimate)
{
    if (Nt <= 2) throw "Nt <= 2";
    if (Nt % 2 == 1) throw "Nt is not a multiple of 2";
    fft_pulse = new fftw_complex[Nt];
    fft_tmp = new fftw_complex[Nt];
    MakeFFTWplans();
    initialised = false;
}

void SimplePulse::CleanUp()
{
    delete[] fft_pulse;
    delete[] fft_tmp;
    fftw_destroy_plan(forward_plan);
    fftw_destroy_plan(backward_plan);
}

SimplePulse::~SimplePulse()
{
    CleanUp();
}

void SimplePulse::MakeFFTWplans()
{
    if (fft_estimate)
    {
	forward_plan = fftw_create_plan(Nt, FFTW_BACKWARD, FFTW_ESTIMATE);
	backward_plan = fftw_create_plan(Nt, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    else
    {
	forward_plan = fftw_create_plan(Nt, FFTW_BACKWARD, FFTW_MEASURE);
	backward_plan = fftw_create_plan(Nt, FFTW_FORWARD, FFTW_MEASURE);
    }
    if (! forward_plan || ! backward_plan)
	throw "ERROR initialising FFTW plans";
}

SimplePulse::SimplePulse(const SimplePulse& p) :
    AbstractPulse(p.Nt, p.t_step, p.is_in_time_domain)
{
    int i;
    fft_pulse = new fftw_complex[Nt];
    fft_tmp = new fftw_complex[Nt];
    fft_estimate = p.fft_estimate;
    MakeFFTWplans();
    if (p.initialised)
    {
	for (i=0; i<Nt; i++)
	{
	    fft_pulse[i].re = p.fft_pulse[i].re;
	    fft_pulse[i].im = p.fft_pulse[i].im;
	}
    }
    initialised = p.initialised;
}


const SimplePulse& SimplePulse::operator= (const SimplePulse& p)
{
    int i;
    if (Nt != p.Nt)
    {
	CleanUp();
	Nt = p.Nt;
	fft_pulse = new fftw_complex[Nt];
	fft_tmp = new fftw_complex[Nt];
	fft_estimate = p.fft_estimate;
	MakeFFTWplans();
    }
    else if (fft_estimate != p.fft_estimate)
    {
	fftw_destroy_plan(forward_plan);
	fftw_destroy_plan(backward_plan);
	fft_estimate = p.fft_estimate;
	MakeFFTWplans();
    }
    t_step = p.t_step;
    if (p.initialised)
    {
	for (i=0; i<Nt; i++)
	{
	    fft_pulse[i].re = p.fft_pulse[i].re;
	    fft_pulse[i].im = p.fft_pulse[i].im;
	}
    }
    initialised = p.initialised;
    return *this;
}


void SimplePulse::ToTimeDomain()
{
    fftw_complex* helper;
    if (initialised)
    {
	if (is_in_time_domain) return;
	fftw_one(backward_plan, fft_pulse, fft_tmp);
	helper = fft_pulse;
	fft_pulse = fft_tmp;
	fft_tmp = helper;
	// FFTW makes a non-normalised Fourier transform -- we have to divide
	// the result by Nt in order to satisfy FFT(IFFT)=1
	for (int i=0; i<Nt; i++)
	{
	    fft_pulse[i].re /= Nt;
	    fft_pulse[i].im /= Nt;
	}
    }
    is_in_time_domain = true;
}

void SimplePulse::ToFrequencyDomain()
{
    fftw_complex* helper;
    if (initialised)
    {
	if (! is_in_time_domain) return;
	fftw_one(forward_plan, fft_pulse, fft_tmp);
	helper = fft_pulse;
	fft_pulse = fft_tmp;
	fft_tmp = helper;
    }
    is_in_time_domain = false;
}

void SimplePulse::DestructiveDomainSwitch(bool to_time_domain)
{
    is_in_time_domain = to_time_domain;
}

void SimplePulse::Get(Complex* array) const
{
    int i;
    if (! initialised) throw "not initialised yet";
    for (i=0; i<Nt/2; i++)
	array[i+Nt/2] = Complex(fft_pulse[i].re, fft_pulse[i].im);
    for (i=Nt/2; i<Nt; i++)
	array[i-Nt/2] = Complex(fft_pulse[i].re, fft_pulse[i].im);
}

void SimplePulse::Set(const Complex* array)
{
    int i;
    Complex z;
    for (i=0; i<Nt/2; i++)
    {
	z = array[i+Nt/2];
	fft_pulse[i].re = real(z);
	fft_pulse[i].im = imag(z);
    }
    for (i=Nt/2; i<Nt; i++)
    {
	z = array[i-Nt/2];
	fft_pulse[i].re = real(z);
	fft_pulse[i].im = imag(z);
    }
    initialised = true;
}

void SimplePulse::GetIntensity(Real* intensity) const
{
    int i;
    Real factor;
    for (i=0; i<Nt/2; i++)
    {
	intensity[i+Nt/2] = square(fft_pulse[i].re) + square(fft_pulse[i].im);
    }
    for (i=Nt/2; i<Nt; i++)
    {
	intensity[i-Nt/2] = square(fft_pulse[i].re) + square(fft_pulse[i].im);
    }
    if (! is_in_time_domain)
    {
	factor = square(t_step);
	for (i=0; i<Nt; i++) intensity[i] *= factor;
    }
}

Real SimplePulse::PeakIntensity() const
{
    int i;
    Real I, Imax=0;
    for (i=0; i<Nt; i++)
    {
	I = square(fft_pulse[i].re) + square(fft_pulse[i].im);
	if (Imax < I) Imax = I;
    }
    if (! is_in_time_domain) Imax *= square(t_step);
    return Imax;
}

void SimplePulse::GetPhase(Real* phase) const
{
    int i;
    Real y;
    for (i=0; i<Nt/2; i++)
    {
	y = square(fft_pulse[i].re) + square(fft_pulse[i].im);
	if (y < 1e-15) phase[i+Nt/2] = 0.0;
	else phase[i+Nt/2] = atan2(fft_pulse[i].im, fft_pulse[i].re);
    }
    for (i=Nt/2; i<Nt; i++)
    {
	y = square(fft_pulse[i].re) + square(fft_pulse[i].im);
	if (y < 1e-15) phase[i-Nt/2] = 0.0;
	else phase[i-Nt/2] = atan2(fft_pulse[i].im, fft_pulse[i].re);
    }
    SmoothenPhase(phase, Nt);
}


/*
void SimplePulse::ACF(Real omega0, Real* acf)
{
    int i;
    Real t, arg;
    fftw_real* A = new fftw_real[Nt];
    fftw_complex* B = new fftw_complex[Nt];
    fftw_complex* C = new fftw_complex[Nt];
    fftw_real* G = new fftw_real[Nt];
    bool was_in_frequency_domain = IsInFrequencyDomain();

    ToTimeDomain();

//     Calculate the 'A' component of the autocorrelation:
//     A(t) = 2*A1(0) + 4*A1(t), where A1(t) = InverseFFT(|FFT(|pulse(t)|^2)|^2).
//     we write temporarily |pulse(t)|^2 into C -- we will need it more than
//     once; we also use B for temporary storage
    for (i=0; i<Nt; i++)
    {
	C[i].re = square(fft_pulse[i].re)+ square(fft_pulse[i].im);
	C[i].im = 0;
    }
    fftw_one(forward_plan, C, fft_tmp);
    for (i=0; i<Nt; i++)
    {
	fft_tmp[i].re = square(fft_tmp[i].re) + square(fft_tmp[i].im);
	fft_tmp[i].im = 0;
    }
    fftw_one(backward_plan, fft_tmp, B);
    for (i=0; i<Nt; i++) A[i] = 4.0*B[i].re + 2.0*B[0].re;

//     Calculate the 'B' component of the autocorrelation:
//     2*InverseFFT(Re[Conjugate(FFT(Conjugate(pulse)*|pulse|^2)) * 
//     FFT(Conjugate(pulse))])
    // we will use B and C for temporary storage;
    // at the moment the real part of C contains |pulse(t)|^2
    for (i=0; i<Nt; i++)
    {
	fft_tmp[i].re = fft_pulse[i].re * C[i].re;
	fft_tmp[i].im = - fft_pulse[i].im * C[i].re;
    }
    fftw_one(forward_plan, fft_tmp, B);
    for (i=0; i<Nt; i++)
    {
	fft_tmp[i].re = fft_pulse[i].re;
	fft_tmp[i].im = - fft_pulse[i].im;
    }
    fftw_one(forward_plan, fft_tmp, C);
    for (i=0; i<Nt; i++)
    {
	fft_tmp[i].re = 2.0 * (B[i].re * C[i].re + 
			       B[i].im * C[i].im);
	fft_tmp[i].im = 0;
    }
    fftw_one(backward_plan, fft_tmp, B);

//     Calculate the 'C' component of the autocorrelation:
//     InverseFFT(|FFT((Conjugate(pulse))^2)|^2)
    for (i=0; i<Nt; i++)
    {
	C[i].re = square(fft_pulse[i].re) - square(fft_pulse[i].im);
	C[i].im = - 2.0 * fft_pulse[i].re * fft_pulse[i].im;
    }
    fftw_one(forward_plan, C, fft_tmp);
    for (i=0; i<Nt; i++)
    {
	fft_tmp[i].re = square(fft_tmp[i].re) + square(fft_tmp[i].im);
	fft_tmp[i].im = 0;
    }
    fftw_one(backward_plan, fft_tmp, C);

    // collect it all into 'G'
    for (i=0; i<Nt; i++)
    {
	if (i < Nt/2) t = t_step * i;
	else t = - t_step * (Nt - i);
	arg = omega0 * t;
	G[i] = A[i] + 
	    4.0 * (B[i].re*cos(arg) - B[i].im*sin(arg)) +
	    2.0 * (C[i].re*cos(2.0*arg) - C[i].im*sin(2.0*arg));
    }

    // normalise the autocorrelation so, that G(0) = 8
    arg = 8.0 / G[0];
    for (i=0; i<Nt; i++)
    {
	G[i] *= arg;
//     A[i] *= arg;
//     B[i].re *= arg;
//     B[i].im *= arg;
//     C[i].re *= arg;
//     C[i].im *= arg;
    }

    // convert from the FFT format to the contiguous format
    for (i=0; i<Nt/2; i++) acf[i+Nt/2] = G[i];
    for (i=Nt/2; i<Nt; i++) acf[i-Nt/2] = G[i];

    if (was_in_frequency_domain) ToFrequencyDomain();

    delete[] A;
    delete[] B;
    delete[] C;
    delete[] G;
}
*/

//-------------------------------------------------------------------------


Pulse::Pulse(int Nt, Real t_step, Real t_min, bool in_time_domain,
	     bool fft_estimate) :
    SimplePulse(Nt, t_step, in_time_domain, fft_estimate), t_min(t_min)
{}


Pulse::Pulse(const Pulse& p) : 
    SimplePulse(p.Nt, p.t_step, p.is_in_time_domain, p.fft_estimate)
{
    SimplePulse::operator=(p);
    t_min = p.t_min;
}

const Pulse& Pulse::operator= (const Pulse& p)
{
    SimplePulse::operator=(p);
    t_min = p.t_min;
    return *this;
}

Pulse::~Pulse() {}

Real Pulse::StartingTime() const
{
    return t_min;
}

Real Pulse::Compare(const SimplePulse& p) const
{
    Complex z;
    Complex* Z = new Complex[Nt];
    int i;
    Real norm1, norm2;
    Real* intensity = new Real[Nt];

    if (this == &p) return 0.0;
    if (Nt != p.NumberOfPoints())
	throw "can't compare pulses of different lengths";
    if (is_in_time_domain && ! p.IsInTimeDomain() ||
	! is_in_time_domain && p.IsInTimeDomain() ) // XOR
	throw "the pulses being compared must be in the same domain";
    p.Get(Z);
    for (i=0; i<Nt/2; i++)
    {
	Z[i+Nt/2] *= Complex(fft_pulse[i].re, -fft_pulse[i].im);
    }
    for (i=Nt/2; i<Nt; i++)
    {
	Z[i-Nt/2] *= Complex(fft_pulse[i].re, -fft_pulse[i].im);
    }
    z = IntegrateArray(Nt, Z, 1.0);
    GetIntensity(intensity);
    norm1 = IntegrateArray(Nt, intensity, 1.0);
    p.GetIntensity(intensity);
    norm2 = IntegrateArray(Nt, intensity, 1.0);
    delete[] Z;
    delete[] intensity;
    return 1.0 - abs(z)/sqrt(norm1*norm2);
}

Real Pulse::CompareIntensities(const SimplePulse& p) const
{
    Real S1, S2, S3;
    int i;
    Real* I1 = new Real[Nt];
    Real* I2 = new Real[Nt];
    GetIntensity(I1);
    p.GetIntensity(I2);
    S1 = 0.5 * (fabs(I1[0] - I2[0]) + fabs(I1[Nt-1] - I2[Nt-1]));
    S2 = 0.5 * (fabs(I1[0]) + fabs(I1[Nt-1]));
    S3 = 0.5 * (fabs(I2[0]) + fabs(I2[Nt-1]));
    for (i = 1; i < Nt-1; i++)
    {
	S1 += fabs(I1[i] - I2[i]);
	S2 += fabs(I1[i]);
	S3 += fabs(I2[i]);
    }
    delete[] I1;
    delete[] I2;
    if (S2 == 0 && S3 == 0) return 0.0;
    return S1 / (S2 + S3);
}

Real Pulse::MeanSquareWidth () const
{
    int i;
    Real s0 = 0;
    Real s1 = 0;
    Real s2 = 0;
    Real x, dx, helper;
    if (is_in_time_domain) dx = TimeStep();
    else dx = OmegaStep();
    for (i=0; i<Nt/2; i++)
    {
	x = dx*i;
	helper = square(fft_pulse[i].re) + square(fft_pulse[i].im);
	s0 += helper;
	s1 += x * helper;
	s2 += x*x * helper;
    }
    for (i=Nt/2; i<Nt; i++)
    {
	x = - dx * (Nt - i);
	helper = square(fft_pulse[i].re) + square(fft_pulse[i].im);
	s0 += helper;
	s1 += x * helper;
	s2 += x*x * helper;
    }
    return sqrt(s2/s0 - square(s1/s0));
}

Real Pulse::FWHM() const
{
    int i, i0, i1, i2;
    Real dx, x1, x2, Imax, dI;
    Real* intensity = new Real[Nt];

    // calculate the intensity, find the maximal intensity (Imax)
    GetIntensity(intensity);
    Imax = 0.0;
    i0 = 0;
    for (i=0; i<Nt; i++)
    {
	if (Imax < intensity[i])
	{
	    Imax = intensity[i];
	    i0 = i;
	}
    }
    // find the points, at which the intensity drops below Imax/2
    for (i1=i0-1; i1>=0 && intensity[i1] >= Imax/2.0; i1--) {}
    if (i1 < 0) i1 = 0;
    for (i2=i0+1; i2<Nt && intensity[i2] >= Imax/2.0; i2++) {}
    if (i2 >= Nt) i2 = Nt - 1;
    // use the linear interpolation to find the values of x, at which
    // the intensity is equal to Imax/2
    if (is_in_time_domain) dx = TimeStep();
    else dx = OmegaStep();
    dI = intensity[i1+1] - intensity[i1];
    if (dI == 0.0) x1 = 0.0;
    else x1 = dx * (i1 + (Imax/2.0-intensity[i1]) / dI);
    dI = intensity[i2-1] - intensity[i2];
    if (dI == 0.0) x2 = 0.0;
    else x2 = dx * (i2 - (Imax/2.0-intensity[i2]) / dI);
    // deallocate 'intensity' and return the result
    delete[] intensity;
    return (x2 - x1);
}

Real Pulse::Energy() const
{
    Real dx, energy;
    Real* intensity = new Real[Nt];
    GetIntensity(intensity);
    if (is_in_time_domain) dx = TimeStep();
    else dx = OmegaStep();
    energy = IntegrateArray(Nt, intensity, dx);
    delete[] intensity;
    return energy;
}

bool Pulse::isFinite() const
{
    int i;
    for (i = 0; i < Nt; i++)
    {
//     if (! isfinite(fft_pulse[i].re) ||
// 	! isfinite(fft_pulse[i].re)) return false;
	if (isinf(fft_pulse[i].re) || isinf(fft_pulse[i].im) ||
	    isnan(fft_pulse[i].re) || isnan(fft_pulse[i].im)) return false;
    }
    return true;
}


Real Pulse::CentreOfGravity() const
{
    int i;
    Real s0 = 0.0;
    Real s1 = 0.0;
    Real x, y;
    // find the centre of the pulse
    for (i=0; i<Nt/2; i++)
    {
	x = (i+Nt/2)*t_step;
	y = square(fft_pulse[i].re) + square(fft_pulse[i].im);
	s0 += y;
	s1 += y * x;
    }
    for (i=Nt/2; i<Nt; i++)
    {
	x = (i-Nt/2)*t_step;
	y = square(fft_pulse[i].re) + square(fft_pulse[i].im);
	s0 += y;
	s1 += y * x;
    }
    return t_min + s1/s0;
}

Real Pulse::PeakPosition() const
{
    int i, i0=0;
    Real t0, I, I_max=0;
    for (i=0; i<Nt; i++)
    {
	I = square(fft_pulse[i].re) + square(fft_pulse[i].im);
	if (I_max < I)
	{
	    I_max = I;
	    i0 = i;
	}
    }
    if (i0<Nt/2) t0 = (i0+Nt/2)*t_step;
    else t0 = (i0-Nt/2)*t_step;
    return t_min + t0;
}

Real Pulse::ShiftToCentre(Real omega0)
{
    int i, j, k;
    Complex phase_correction;
    Complex* Z = new Complex[Nt];
    Real shift;
    bool was_in_time_domain = is_in_time_domain;

    ToTimeDomain();
    // convert from the FFT format to the contiguous format
    Get(Z);
    // determine the amount of shift -> j
    j = (int) rint((Nt-1)/Real(2) - (PeakPosition()-t_min)/t_step);
    shift = j*t_step;
    phase_correction = polar(Real(1),omega0*shift);
    t_min -= shift;
    if (j>0) // shift the pulse to the right
    {
	for (i=Nt-1, k=i-j; k>=0; i--, k--) Z[i] = Z[k] * phase_correction;
	// pad with zeros
	for (i=0; i<j; i++) Z[i] = 0.0;
    }
    else if (j<0) // shift the pulse to the left
    {
	j = -j;
	for (i=0, k=j; k<Nt; i++, k++) Z[i] = Z[k] * phase_correction;
	// pad with zeros
	for (i = Nt-j; i<Nt; i++) Z[i] = 0.0;
    }
    Set(Z);
    if (! was_in_time_domain) ToFrequencyDomain();
    delete[] Z;
    return shift;
}

void Pulse::Normalise(Real I_max)
{
    int i;
    Real I0, factor;
    Real* intensity = new Real[Nt];

    // calculate the intensity, find the maximal intensity -> I0
    GetIntensity(intensity);
    I0 = 0.0;
    for (i=0; i<Nt; i++)
	if (I0 < intensity[i]) I0 = intensity[i];
    if (I0 > 0.0)
    { // normalise
	factor = sqrt(I_max/I0);
	for (i=0; i<Nt; i++)
	{
	    fft_pulse[i].re *= factor;
	    fft_pulse[i].im *= factor;
	}
    }
    delete[] intensity;
}
