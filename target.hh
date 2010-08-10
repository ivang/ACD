// optimisation targets

#ifndef __TARGET_HH
#define __TARGET_HH

#include "common.hh"
#include "parameters.hh"
#include "coating.hh"
#include "pulse.hh"

class AbstractTarget
{
public:
    AbstractTarget(int number_of_frequencies, Real omega_min, Real omega_max,
		   int number_of_bounces,
		   Real minimal_layer_thickness=0.0,
		   Real individualism=0.0);
   /* The discretisation with respect to the frequency is fixed and
      MAY NOT BE CHANGED DURING THE LIFETIME OF THE OBJECT
      (this is the price for having 'MeritFunction' and
      'MeritGradient' implemented at this level of abstraction).
      The parameter 'individualism' is a real number in the range [0,1],
      which is used for the simultaneous optimisation of a set of
      complementary mirrors: if it is equal to zero, then only the
      average properties of the set matter, if it is equal to 1, then
      the average properties are not considered at all--the sum of
      the figures of merit is minimised. */

    virtual ~AbstractTarget();

    virtual Real NonReflectivityMerit(Coating& coating,
				      int perturbed_layer=-1,
				      Real thickness=0.0);
    virtual Real NonReflectivityMerit(vector<Coating>& coatings,
				      int perturbed_coating=-1,
				      int perturbed_layer=-1,
				      Real thickness=0.0);
    /* The figure of merit is not always completely determined by
       the reflectivity; sometimes other coating parameters play an
       important role. An example of such a factor is a constraint to
       possible layer thicknesses. The implementations of these two
       functions in this class return some "punishment" value, if one
       of the layers is too thin (this value is automatically
       added to the figure of merit). Children of 'AbstractTarget' 
       may extend or substitute this behaviour. */

    virtual Real ThinLayerPunishment(Real thickness);
    // defines the punishment strategy

    // evaluate the merit function of a coating
    virtual Real MeritFunction(Coating& coating);

    // evaluate the merit function of a coating, assuming a different
    // thickness for one of the layers (we can do it fast)
    virtual Real MeritFunction(Coating& coating, int layer,
			       Real thickness);

    // evaluate the merit function of a set of complementary mirrors
    virtual Real MeritFunction(vector<Coating>&);

    // evaluate the gradient of the merit function
    virtual void MeritGradient(Coating& coating, Real* gradient);
    virtual void MeritGradient(vector<Coating>&, Real* gradient);
    Real small_thickness_variation; // the thickness variation used
				    // to calculate the gradient
protected:
    int number_of_frequencies;
    Real omega_min, omega_max;
    int number_of_bounces;
    Real minimal_layer_thickness;
    Real individualism;
    Real *frequencies;
    Real *wavelengths;
    Complex *reflectivity;
    Real* reflectance;
    Real* phase_shift;

    virtual Real Merit(Coating& coating, Complex* reflectivity=NULL, 
	    Real* _reflectance=NULL, Real* _phase_shift=NULL) = 0;
    /* this function is the core of the abstraction: different targets
       calculate the figure of merit in different ways; the
       number of frequencies from 'number_of_frequencies' is assumed to
       be the size of the array 'reflectivity' */

    void CheckCoating(const Coating& coating);
    void CheckCoatings(const vector<Coating>& coatings);
};

//-------------------------------------------------------------------------

class GDDTarget: public AbstractTarget
{
 public:
    GDDTarget(int number_of_frequencies, Real omega_min, Real omega_max,
	      AbstractDispersion& target_dispersion,
	      AbstractParameter& target_reflectance,
	      AbstractParameter& GD_tolerance,
	      AbstractParameter& GD_reserve,
	      AbstractParameter& GDD_tolerance,
	      AbstractParameter& GDD_reserve,
	      Real GD_contribution,
	      AbstractParameter& reflectance_tolerance,
	      AbstractParameter& reflectance_reserve,
	      int number_of_bounces,
	      Real minimal_layer_thickness=0.0,
	      Real individualism=0.0,
	      int merit_power=2,
	      bool adaptive_dispersion=false);
    ~GDDTarget();

protected:
    Real* target_GD;
    Real* adapted_target_GD;
    Real* target_GDD;
    Real* adapted_target_GDD;
    Real* target_reflectance;
    Real* GD_tolerance;
    Real* GD_reserve;
    Real* GDD_tolerance;
    Real* GDD_reserve;
    Real* reflectance_tolerance;
    Real* reflectance_reserve;
    Real GD_contribution;
    int merit_power;
    bool adaptive_dispersion;
    Real* GD;
    Real* GDD;
    Real Merit(Coating& coating, Complex* reflectivity=NULL,
	    Real* _reflectance=NULL, Real* _phase_shift=NULL);
};

//-------------------------------------------------------------------------

class PulseTarget: public AbstractTarget
{
 public:
    PulseTarget(int number_of_frequencies, Real omega_min, Real omega_max,
		AbstractParameter& probe_pulse_spectrum,
		AbstractParameter& probe_pulse_phase,
		AbstractDispersion& target_dispersion,
		Real suggested_time_resolution, // in femtoseconds
		int number_of_bounces,
		Real minimal_layer_thickness=0.0,
		Real individualism=0.0,
		bool adaptive_dispersion=false);
    ~PulseTarget();

protected:
    bool adaptive_dispersion;
    Real* target_phase_shift; // used to handle adaptive dispersion
    Real* target_GDD; // used to handle adaptive dispersion
    Real* GDD; // used to handle adaptive dispersion
    int N_pulse; // the number of points for the pulse
    Complex* probe_pulse; // prechirped with the dispersion, which the mirror
                          // has to compensate for
    Real* AD_weights; // used to handle adaptive dispersion
    Complex* reflected_pulse;
    Pulse* pulse;
    Real Merit(Coating& coating, Complex* reflectivity=NULL,
	    Real* _reflectance=NULL, Real* _phase_shift=NULL);

    Real AD_Factor(const Real* GDD) const;
    /* an analog of 'AdaptiveDispersionFactor', which uses the given weights
       and doesn't know anything about GDD_reserve */

};

#endif

