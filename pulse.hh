#ifndef __PULSE_HH
#define __PULSE_HH

#include "common.hh"

#ifdef SINGLE_PRECISION
#include "sfftw.h"
#else
#include "fftw.h"
#endif


// The pulse is supposed to be specified by its envelope

class AbstractPulse
{
protected:
    int Nt;
    Real t_step;
    bool is_in_time_domain;
public:
    AbstractPulse(int Nt, Real t_step, bool in_time_domain=true);
    virtual ~AbstractPulse();
    int NumberOfPoints() const;
    Real TimeStep() const;
    Real OmegaStep() const;
    virtual void ToTimeDomain() = 0;
    virtual void ToFrequencyDomain() = 0;
    inline bool IsInTimeDomain() const;
    inline bool IsInFrequencyDomain() const;
    virtual void Get(Complex* array) const = 0; // the array must have Nt elements
    virtual void Set(const Complex* array) = 0;
    // the array must have Nt elements
};

class SimplePulse : public AbstractPulse
{
private:
    fftw_plan forward_plan, backward_plan;
    void MakeFFTWplans();
    bool initialised;
    void CleanUp();
protected:
    fftw_complex *fft_pulse; // the pulse in the current domain in the
    // FFT format
    fftw_complex *fft_tmp; // for temporary storage
    bool fft_estimate;
public:
    SimplePulse(int Nt, Real t_step, bool in_time_domain=true,
		bool fft_estimate=true);
    /* if 'fft_estimate' is true, the FFT initialisation performed by
       the constructor is fast, but the consequent Fourier transforms
       will be somewhat slow; if 'fft_estimate' is false, the constructor
       may take a few seconds to initialise FFT, but the the FFT evaluations
       are very fast */

    SimplePulse(const SimplePulse&);
    ~SimplePulse();
    const SimplePulse& operator= (const SimplePulse&);
    void ToTimeDomain();
    void ToFrequencyDomain();

    void DestructiveDomainSwitch(bool to_time_domain);
    // makes no FFT, just changes the current domain

    void Get(Complex*) const;
    void Set(const Complex*);
    void GetIntensity(Real*) const; // writes |A|^2 to the given array
    Real PeakIntensity() const;
    void GetPhase(Real*) const;

//    void ACF(Real omega0, Real* acf);
    /* calculates the second-order autocorrelation function (IACF) in the
       time domain; the array 'acf' must have NumberOfPoints() elements,
       it will be filled with acf[i] = acf(dt*(i-N/2)) */
};

class Pulse : public SimplePulse
{
protected:
    Real t_min;
public:
    Pulse(int Nt, Real t_step, Real t_min=0.0,
	  bool in_time_domain=true, bool fft_estimate=true);
    Pulse(const Pulse&);
    const Pulse& operator= (const Pulse&);
//   Pulse(const SimplePulse&, Real tmin);
//   const Pulse& operator= (const SimplePulse&, Real tmin);
    ~Pulse();

    Real StartingTime() const; // returns t_min -- the time corresponding
    // to the very first element

    Real Compare(const SimplePulse&) const;
    /* compares two objects of the class; returns a value between 0 and 1;
       0 indicates that the two pulses are identical, 1 indicates
       that the pulses are absolutely different; THIS FUNCTION IGNORES
       THE CARRIER-ENVELOPE PHASE */

    Real CompareIntensities(const SimplePulse&) const;
    /* like 'Compare', but compares only the intensities of the pulses,
       if they are in the time domain, or their power spectra, if they
       are in the frequency domain. */

    Real MeanSquareWidth () const;
    /* Returns the mean square width of the pulse intensity, if the pulse
       is in the time domain, or the RMS of the power spectrum, if the
       pulse is in the frequency domain */

    bool isFinite() const; // checks that the complex amplitude is finite
    Real FWHM() const; // returns the FWHM of the INTENSITY
    Real Energy() const; // integrates the intensity
    Real CentreOfGravity() const;
    Real PeakPosition() const;
    Real ShiftToCentre(Real omega0); // returns the shift
    void Normalise(Real I_max=1.0); // makes max[|A|^2]=I_max
    //  bool isTooCloseToBoundary() const;
};

#endif
