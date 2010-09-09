#include "target.hh"
#include <cmath>

//-------------------------------------------------------------------------

AbstractTarget::AbstractTarget(int number_of_frequencies,
			       Real omega_min, Real omega_max,
			       int number_of_bounces,
			       Real minimal_layer_thickness,
			       Real individualism) :
    number_of_frequencies(number_of_frequencies),
    omega_min(omega_min), omega_max(omega_max),
    number_of_bounces(number_of_bounces),
    minimal_layer_thickness(minimal_layer_thickness),
    individualism(individualism)
{
    if (number_of_frequencies < 2)
    {
	throw "ERROR in AbstractTarget::AbstractTarget:\n\
number_of_frequencies < 2";
    }
    small_thickness_variation = 0.1; // the default value
    frequencies = new Real[number_of_frequencies];
    wavelengths = new Real[number_of_frequencies];
    reflectivity = new Complex[number_of_frequencies];
    reflectance = new Real[number_of_frequencies];
    phase_shift = new Real[number_of_frequencies];
    for (int i=0; i<number_of_frequencies; i++)
    {
	frequencies[i] = omega_min +
	    (omega_max - omega_min) * i/Real(number_of_frequencies-1);
	wavelengths[i] = 2.0*M_PI*SPEED_OF_LIGHT / frequencies[i];
    }
}


AbstractTarget::~AbstractTarget()
{
    delete[] frequencies;
    delete[] wavelengths;
    delete[] reflectivity;
    delete[] reflectance;
    delete[] phase_shift;
}

//-------------------------------------------------------------------------

void AbstractTarget::CheckCoating(const Coating& coating)
{
    Coating::Parameters parameters;
    int N = coating.NumberOfFrequencies();
    coating.GetParameters(parameters);
    if (N != number_of_frequencies ||
	parameters.omega_min != omega_min ||
	parameters.omega_max != omega_max)
    {
	throw "ERROR in AbstractTarget::CheckCoating:\n\
a different discretisation is expected";
    }
}

//-------------------------------------------------------------------------

void AbstractTarget::CheckCoatings(const vector<Coating>& coatings)
{
    int N;
    Coating::Parameters parameters;
    int n = coatings.size();
    if (n <= 0) throw "ERROR in AbstractTarget::CheckCoatings: \
the vector is empty";
    for (int i=0; i<n; i++)
    {
	N = coatings[i].NumberOfFrequencies();
	coatings[i].GetParameters(parameters);
	if (N != number_of_frequencies ||
	    parameters.omega_min != omega_min ||
	    parameters.omega_max != omega_max)
	    throw "ERROR in AbstractTarget::CheckCoatings: \
a different discretisation is expected";
    }
}

//-------------------------------------------------------------------------

Real AbstractTarget::MeritFunction(Coating& coating)
{
    Real merit;
    CheckCoating(coating);
    coating.Reflectivity(reflectivity);
    merit = Merit(coating) + NonReflectivityMerit(coating);
    if (isinf(merit) || isnan(merit))
    {
	throw "ERROR in AbstractTarget::MeritFunction:\n\
the figure of merit is not a finite number";
    }
    return merit;
}

//-------------------------------------------------------------------------

Real AbstractTarget::MeritFunction(Coating& coating, int layer,
				   Real thickness)
{
    Real merit;
    CheckCoating(coating);
    coating.TryLayerThickness(layer, thickness, reflectivity);
    merit = Merit(coating, reflectivity);
    merit += NonReflectivityMerit(coating, layer, thickness);
    if (isinf(merit) || isnan(merit))
    {
	throw "ERROR in AbstractTarget::MeritFunction:\n\
the figure of merit is not a finite number";
    }
    return merit;
}

//-------------------------------------------------------------------------

Real AbstractTarget::MeritFunction(vector<Coating>& coatings)
{   // the merit function of a set of complementary mirrors
    int i, k;
    Real sum_of_merits, merit_of_the_set, merit;
    Real *net_reflectance, *net_phase_shift;
    int n = coatings.size();
    CheckCoatings(coatings);
    // analyse all the mirrors one by one and accumulate the
    // information about the set as a whole
    net_reflectance = new Real[number_of_frequencies];
    net_phase_shift = new Real[number_of_frequencies];
    for (i=0; i<number_of_frequencies; i++)
    {
	net_reflectance[i] = 1;
	net_phase_shift[i] = 0;
    }
    sum_of_merits = 0.0;
    for (k=0; k<n; k++)
    {
	coatings[k].Reflectivity(reflectivity);
	Reflectance(number_of_frequencies, reflectivity, reflectance);
	PhaseShift(number_of_frequencies, reflectivity, phase_shift);
	for (i=0; i<number_of_frequencies; i++)
	{
	    net_reflectance[i] *= reflectance[i];
	    net_phase_shift[i] += phase_shift[i];
	}
	if (individualism != Real(0))
	    sum_of_merits += Merit(coatings[k], reflectivity);
    }
    // analyse the set of mirrors as a whole
    for (i=0; i<number_of_frequencies; i++)
    {
	net_reflectance[i] = pow(net_reflectance[i], Real(1)/n);
	net_phase_shift[i] /= n;
    }
    merit_of_the_set = Merit(coatings[1], reflectivity=NULL,
	    net_reflectance, net_phase_shift);
    merit = individualism*sum_of_merits/n + 
	(Real(1)-individualism)*merit_of_the_set;
    merit += NonReflectivityMerit(coatings);
    delete[] net_reflectance;
    delete[] net_phase_shift;
    if (isinf(merit) || isnan(merit))
    {
	throw "ERROR in AbstractTarget::MeritFunction:\n\
the figure of merit is not a finite number";
    }
    return merit;
}

//-------------------------------------------------------------------------

void AbstractTarget::MeritGradient(Coating& coating, Real* gradient)
{
    int i, n;
    Real d, f1, f2;
    Real* thicknesses;

    CheckCoating(coating);
    n = coating.NumberOfVariableLayers();
    if (n<=0) return;
    thicknesses = new Real[n];
    coating.GetVariableThicknesses(thicknesses);
    for (i=0; i<n; i++)
    {
	d = thicknesses[i] + small_thickness_variation;
	coating.TryLayerThickness(i, d, reflectivity);
	f2 = Merit(coating, reflectivity) +
	    NonReflectivityMerit(coating, i, d);
	d = thicknesses[i] - small_thickness_variation;
	coating.TryLayerThickness(i, d, reflectivity);
	f1 = Merit(coating, reflectivity) +
	    NonReflectivityMerit(coating, i, d);
	gradient[i] = (f2 - f1) / (2.0*small_thickness_variation);
    }
    delete[] thicknesses;
}

//-------------------------------------------------------------------------

void AbstractTarget::MeritGradient(vector<Coating>& coatings, Real* gradient)
{
    int i, j, m, n, coating_index, k, number_of_coatings, N;
    Real d, f1, f2, merit;
    Real *thicknesses=NULL, *merits=NULL;
    Real *net_reflectance, *net_phase_shift;
    Complex** reflectivities;
    Real sum_of_merits, merit_of_the_set;

    CheckCoatings(coatings);
    N = number_of_frequencies;
    number_of_coatings = coatings.size();
    // calculate the reflectivity and the figure of merit for
    // each of the coatings
    merits = new Real[number_of_coatings];
    net_reflectance = new Real[N];
    net_phase_shift = new Real[N];
    reflectivities = new Complex*[number_of_coatings];
    sum_of_merits = 0;
    for (coating_index=0; coating_index<number_of_coatings; coating_index++)
    {
	reflectivities[coating_index] = new Complex[N];
	coatings[coating_index].Reflectivity(reflectivity);
	for (i=0; i<N; i++) reflectivities[coating_index][i] = reflectivity[i];
	merits[coating_index] = Merit(coatings[coating_index], reflectivity);
	sum_of_merits += merits[coating_index];
    }
    // calculate the gradient
    m = 0; // the index of the element in the gradient vector
    for (coating_index=0; coating_index<number_of_coatings; coating_index++)
    {
	n = coatings[coating_index].NumberOfVariableLayers();
	if (n>0)
	{
	    thicknesses = new Real[n];
	    coatings[coating_index].GetVariableThicknesses(thicknesses);
	}
	for (i=0; i<n; i++)
	{
	    d = thicknesses[i] + small_thickness_variation;
	    coatings[coating_index].TryLayerThickness(i, d, reflectivity);
	    // calculate the figure of merit for the perturbed coating
	    merit = Merit(coatings[coating_index], reflectivity);
            // calculate the net reflectivity of the set of mirrors
	    for (j=0; j<N; j++)
	    {
		net_reflectance[j] = 1;
		net_phase_shift[j] = 0;
	    }
	    for (k=0; k<number_of_coatings; k++)
	    {
		if (k==coating_index)
		{
		    Reflectance(N, reflectivity, reflectance);
		    PhaseShift(N, reflectivity, phase_shift);
		    for (j=0; j<N; j++)
		    {
			net_reflectance[j] *= reflectance[j];
			net_phase_shift[j] += phase_shift[j];
		    }
		}
		else
		{
		    Reflectance(N, reflectivities[k], reflectance);
		    PhaseShift(N, reflectivities[k], phase_shift);
		    for (j=0; j<N; j++)
		    {
			net_reflectance[j] *= reflectance[j];
			net_phase_shift[j] += phase_shift[j];
		    }
		}
	    }
	    for (j=0; j<N; j++)
	    {
		net_reflectance[j] =
		    pow(net_reflectance[j], Real(1)/number_of_coatings);
		net_phase_shift[j] /= number_of_coatings;
	    }
	    // analyse the system of mirrors as a whole
	    merit_of_the_set = Merit(coatings[1], reflectivity=NULL,
		    net_reflectance, net_phase_shift);
	    // calculate the perturbed merit function
	    f2 = individualism *
		(sum_of_merits-merits[coating_index] + merit) / 
		number_of_coatings + 
		(1.0-individualism)*merit_of_the_set;
	    f2 += NonReflectivityMerit(coatings, coating_index, i, d);
	    // do the same with the other sign of the perturbations
	    d = thicknesses[i] - small_thickness_variation;
	    coatings[coating_index].TryLayerThickness(i, d, reflectivity);
	    // calculate the figure of merit for the perturbed coating
	    merit = Merit(coatings[coating_index], reflectivity);
            // calculate the net reflectivity of the set of mirrors
	    for (j=0; j<N; j++)
	    {
		net_reflectance[j] = 1;
		net_phase_shift[j] = 0;
	    }
	    for (k=0; k<number_of_coatings; k++)
	    {
		if (k==coating_index)
		{
		    Reflectance(N, reflectivity, reflectance);
		    PhaseShift(N, reflectivity, phase_shift);
		    for (j=0; j<N; j++)
		    {
			net_reflectance[j] *= reflectance[j];
			net_phase_shift[j] += phase_shift[j];
		    }
		}
		else
		{
		    Reflectance(N, reflectivities[k], reflectance);
		    PhaseShift(N, reflectivities[k], phase_shift);
		    for (j=0; j<N; j++)
		    {
			net_reflectance[j] *= reflectance[j];
			net_phase_shift[j] += phase_shift[j];
		    }
		}
	    }
	    for (j=0; j<N; j++)
	    {
		net_reflectance[j] =
		    pow(net_reflectance[j], Real(1)/number_of_coatings);
		net_phase_shift[j] /= number_of_coatings;
	    }
	    // analyse the system of mirrors as a whole
	    merit_of_the_set = Merit(coatings[1], reflectivity=NULL, 
		    net_reflectance, net_phase_shift);
	    // calculate the perturbed merit function
	    f1 = individualism *
		(sum_of_merits-merits[coating_index]+ merit) /
		number_of_coatings + 
		(1.0-individualism)*merit_of_the_set;
	    f1 += NonReflectivityMerit(coatings, coating_index, i, d);
	    // calculate the element of the gradient
	    gradient[m++] = (f2 - f1) / (2.0*small_thickness_variation);
	}
	if (n>0) delete[] thicknesses;
    }
    // deallocate the memory
    for (k=0; k<number_of_coatings; k++) delete[] reflectivities[k];
    delete[] reflectivities;
    delete[] merits;
    delete[] net_reflectance;
    delete[] net_phase_shift;
}

//-------------------------------------------------------------------------

Real AbstractTarget::NonReflectivityMerit(Coating& coating,
					  int perturbed_layer,
					  Real thickness)
{
    Real result = 0.0;
    thickness = abs(thickness);
    if (minimal_layer_thickness > 0.0)
    {
	// find the thinnest layer
	Real min_thickness;
	int n = coating.NumberOfVariableLayers();
	if (n<=0) return 0.0;
	Real* thicknesses = new Real[n];
	coating.GetVariableThicknesses(thicknesses);
	if (perturbed_layer>=0 && perturbed_layer<n)
	    thicknesses[perturbed_layer] = thickness;
	min_thickness = thicknesses[0];
	for (int i=1; i<n; i++)
	{
	    if (min_thickness > thicknesses[i])
		min_thickness = thicknesses[i];
	}
	if (min_thickness < minimal_layer_thickness)
	    result = ThinLayerPunishment(min_thickness);
	delete[] thicknesses;
    }
    return result;
}


Real AbstractTarget::NonReflectivityMerit(vector<Coating>& coatings,
					  int perturbed_coating,
					  int perturbed_layer,
					  Real thickness)
{
    Real result = 0.0;
    thickness = abs(thickness);
    if (minimal_layer_thickness > 0.0)
    {
	int k, m;
	int n = coatings.size();
	// find the thinnest layer
	Real min_thickness = 0.0;
	bool first = true;
	for (k=0; k<n; k++)
	{
	    m = coatings[k].NumberOfVariableLayers();
	    if (m<=0) continue;
	    Real* thicknesses = new Real[m];
	    coatings[k].GetVariableThicknesses(thicknesses);
	    if (k==perturbed_coating && 
		perturbed_layer>=0 && perturbed_layer<m)
		thicknesses[perturbed_layer] = thickness;
	    for (int i=0; i<m; i++)
	    {
		if (first)
		{
		    min_thickness = thicknesses[i];
		    first = false;
		}
		else if (min_thickness > thicknesses[i])
		    min_thickness = thicknesses[i];
	    }
	    delete[] thicknesses;
	}
	if (min_thickness < minimal_layer_thickness)
	    result = ThinLayerPunishment(min_thickness);
    }
    return result;
}


Real AbstractTarget::ThinLayerPunishment(Real thickness)
{
    if (thickness==Real(0)) return 0;
    thickness = abs(thickness); // just to be sure
    return 100.0 * square(minimal_layer_thickness / thickness - 1.0);
}

//-------------------------------------------------------------------------

GDDTarget::GDDTarget(int number_of_frequencies,
		     Real omega_min, Real omega_max,
		     AbstractDispersion& target_dispersion,
		     AbstractParameter& _target_reflectance,
		     AbstractParameter& _GD_tolerance,
		     AbstractParameter& _GD_reserve,
		     AbstractParameter& _GDD_tolerance,
		     AbstractParameter& _GDD_reserve,
		     Real GD_contribution,
		     AbstractParameter& _reflectance_tolerance,
		     AbstractParameter& _reflectance_reserve,
		     int number_of_bounces,
		     Real target_EF,
		     Real EF_tolerance,
		     Real EF_reserve,
		     Real dx,
		     Real minimal_layer_thickness,
		     Real individualism,
		     int merit_power,
		     bool adaptive_dispersion) :
    AbstractTarget(number_of_frequencies, omega_min, omega_max,
		   number_of_bounces, minimal_layer_thickness, individualism),
    GD_contribution(GD_contribution),
    merit_power(merit_power),
    adaptive_dispersion(adaptive_dispersion),
    target_EF(target_EF), 
    EF_tolerance(EF_tolerance), EF_reserve(EF_reserve), dx(dx)
{
    int i;
    GD = new Real[number_of_frequencies];
    GDD = new Real[number_of_frequencies];
    target_GD = new Real[number_of_frequencies];
    target_GDD = new Real[number_of_frequencies];
    target_reflectance = new Real[number_of_frequencies];
    GD_tolerance = new Real[number_of_frequencies];
    GD_reserve = new Real[number_of_frequencies];
    GDD_tolerance = new Real[number_of_frequencies];
    GDD_reserve = new Real[number_of_frequencies];
    reflectance_tolerance = new Real[number_of_frequencies];
    reflectance_reserve = new Real[number_of_frequencies];
    adapted_target_GD = new Real[number_of_frequencies];
    adapted_target_GDD = new Real[number_of_frequencies];
    for (i=0; i<number_of_frequencies; i++)
    {
	target_GD[i] = target_dispersion.GD(frequencies[i]);
	target_GDD[i] = target_dispersion.GDD(frequencies[i]);
	target_reflectance[i] = _target_reflectance(wavelengths[i]);
	GD_tolerance[i] = _GD_tolerance(wavelengths[i]);
	GD_reserve[i] = _GD_reserve(wavelengths[i]);
	GDD_tolerance[i] = _GDD_tolerance(wavelengths[i]);
	GDD_reserve[i] = _GDD_reserve(wavelengths[i]);
	reflectance_tolerance[i] = _reflectance_tolerance(wavelengths[i]);
	reflectance_reserve[i] = _reflectance_reserve(wavelengths[i]);
	if (GD_tolerance[i]<=Real(0) || GDD_tolerance[i]<=Real(0) ||
	    reflectance_tolerance[i]<=Real(0))
	{
	    throw "tolerances must be positive";
	}
    }
    if (! adaptive_dispersion)
    {
	for (i=0; i<number_of_frequencies; i++)
	{
	    adapted_target_GD[i] = target_GD[i];
	    adapted_target_GDD[i] = target_GDD[i];
	}
    }
}

//-------------------------------------------------------------------------

GDDTarget::~GDDTarget()
{
    delete[] GD;
    delete[] GDD;
    delete[] target_GD;
    delete[] target_GDD;
    delete[] target_reflectance;
    delete[] GD_tolerance;
    delete[] GD_reserve;
    delete[] GDD_tolerance;
    delete[] GDD_reserve;
    delete[] reflectance_tolerance;
    delete[] reflectance_reserve;
    delete[] adapted_target_GD;
    delete[] adapted_target_GDD;
}

//-------------------------------------------------------------------------

Real GDDTarget::Merit(Coating& coating, Complex* reflectivity,
       Real* _reflectance, Real* _phase_shift)
{
    int i, k, m, N, Nx;
    Real GD_shift, x, reserve, tolerance, omega_step, result;
    Real merit_from_R, merit_from_GD, merit_from_GDD, merit_from_E;
    Real **Z;

    if (!_reflectance) _reflectance = new Real[number_of_frequencies];
    if (!_phase_shift) _phase_shift = new Real[number_of_frequencies];
    
    if (!reflectivity) {
	reflectivity = new Complex[number_of_frequencies];
	coating.Reflectivity(reflectivity);
    }
    
    Reflectance(number_of_frequencies, reflectivity, _reflectance);
    PhaseShift(number_of_frequencies, reflectivity, _phase_shift);

    N = number_of_frequencies;
    omega_step = (frequencies[N-1] - frequencies[0]) / Real(N-1);
    PhaseToGD(N, omega_step, _phase_shift, GD);
    PhaseToGDD(N, omega_step, _phase_shift, GDD);
    
    // maybe we make more than one bounce
    if (number_of_bounces > 1)
    {
	for (i=0; i<N; i++)
	{
	    _reflectance[i] = IntPower(_reflectance[i], number_of_bounces);
	    GD[i] *= number_of_bounces;
	    GDD[i] *= number_of_bounces;
	}
    }

    if (adaptive_dispersion)
    {
	Real factor =
	    AdaptiveDispersionFactor(N, target_GDD, GDD_reserve,
				     GDD_tolerance, GDD);
	for (i=0; i<N; i++)
	{
	    adapted_target_GD[i] = target_GD[i]*factor;
	    adapted_target_GDD[i] = target_GDD[i]*factor;
	}
    }

    // the absolute group delay is of no interest to us; let's
    // calculate the amount of GD, by which we have to shift the
    // calculated GD in order to make it as close to the target GD
    // as possible
    GD_shift = GDShift(N, adapted_target_GD, GD_reserve, GD_tolerance, GD);

    // analyse the properties
    merit_from_R = 0.0;
    merit_from_GD = 0.0;
    merit_from_GDD = 0.0;
    merit_from_E = 0.0;

    // discretisation of the stack thickness
    Nx = int(ceil(coating.StackThickness()/dx));
    Z = new Real*[Nx];
    
    m = 0;

    // calculate the E-field intensity distribution and store it in 'Z'
    coating.EFieldIntensity(N, Nx, dx, Z);

    for (i=0; i<N; i++)
    {
	if (GD_contribution != Real(0))
	{
	    x = abs(GD[i] + GD_shift - adapted_target_GD[i]);
	    reserve = GD_reserve[i];
	    if (x > reserve)
	    {
		tolerance = GD_tolerance[i];
		merit_from_GD += IntPower((x-reserve)/tolerance, merit_power);
	    }
	}
	if (GD_contribution != Real(1))
	{
	    x = abs(GDD[i] - adapted_target_GDD[i]);
	    reserve = GDD_reserve[i];
	    if (x > reserve)
	    {
		tolerance = GDD_tolerance[i];
		merit_from_GDD += IntPower((x-reserve)/tolerance, merit_power);
	    }
	}
	x = abs(_reflectance[i] - target_reflectance[i]);
	reserve = reflectance_reserve[i];
	if (x > reserve)
	{
	    tolerance = reflectance_tolerance[i];
	    merit_from_R += IntPower((x-reserve)/tolerance, merit_power);
	}
	if (target_EF > Real(0))
	{
	    for (k=0; k<Nx; k++)
	    {
		x = Z[k][i] - target_EF;
		reserve = EF_reserve;
		if (x > reserve && x > 0) {
		    tolerance = EF_tolerance;
		    merit_from_E += IntPower((x-EF_reserve)/EF_tolerance, merit_power);
		    m++;
		}
	    }
	}
    }
    
    if (m > 0) merit_from_E /= m / Nx;
    
    result = merit_from_R + merit_from_GD*GD_contribution +
	merit_from_GDD*(Real(1)-GD_contribution) + merit_from_E;
    result = pow(result/(3*N), Real(1)/merit_power);
    
    for (k=0; k<Nx; k++) delete[] Z[k];
    delete[] Z;

    return result;
}

//=========================================================================

PulseTarget::PulseTarget(int number_of_frequencies,
			 Real omega_min, Real omega_max,
			 AbstractParameter& probe_pulse_spectrum,
			 AbstractParameter& probe_pulse_phase,
			 AbstractDispersion& target_dispersion,
			 Real suggested_time_resolution,
			 int number_of_bounces,
			 Real minimal_layer_thickness,
			 Real individualism,
			 bool adaptive_dispersion) :
    AbstractTarget(number_of_frequencies, omega_min, omega_max,
		   number_of_bounces, minimal_layer_thickness, individualism),
    adaptive_dispersion(adaptive_dispersion)
{
    int i, i1, j;
    Real pulse_dt, pulse_domega, A, phi;
    if (adaptive_dispersion)
    {
	target_phase_shift = new Real[number_of_frequencies];
	target_GDD = new Real[number_of_frequencies];
	GDD = new Real[number_of_frequencies];
	AD_weights = new Real[number_of_frequencies];
	for (i=0; i<number_of_frequencies; i++)
	{
	    target_phase_shift[i] = target_dispersion.Phase(frequencies[i]);
	    target_GDD[i] = target_dispersion.GDD(frequencies[i]);
	    AD_weights[i] = probe_pulse_spectrum(wavelengths[i]);
	}
    }
    else
    {
	target_phase_shift = NULL;
	target_GDD = NULL;
	GDD = NULL;
	AD_weights = NULL;
    }
    // determine the number of points for the probe pulse
    N_pulse = IntPower(2, (int) ceil(log(Real(number_of_frequencies))/
				     log(Real(2))));
    pulse_domega = (omega_max - omega_min) / (number_of_frequencies-1);
    pulse_dt = 2.0*M_PI / (N_pulse*pulse_domega);
    if (suggested_time_resolution <= Real(0))
    {
	throw "ERROR in PulseTarget::PulseTarget: \
suggested_time_resolution <= 0";
    }
    while (pulse_dt > suggested_time_resolution*1.1)
    {
	N_pulse *= 2;
	pulse_dt /= 2.0;
    }
#ifdef DEBUG
    cout << "using " << N_pulse << " points to simulate pulse reflection"
	 << endl;
#endif
    probe_pulse = new Complex[N_pulse];
    reflected_pulse = new Complex[N_pulse];
    pulse = new Pulse(N_pulse, pulse_dt, -pulse_dt*(N_pulse/2), false);
    // calculate the probe pulse (do not prechirp)
    if (N_pulse < number_of_frequencies)
	throw "ERROR in PulseTarget::PulseTarget: \
N_pulse < number_of_frequencies";
    i1 = (N_pulse-number_of_frequencies)/2;
    for (i=0; i<i1; i++) probe_pulse[i] = 0;
    for (i=i1,j=0; j<number_of_frequencies; j++,i++)
    {
	A = sqrt(probe_pulse_spectrum(wavelengths[j]));
	A /= frequencies[j]; // I(omega)domega = I(lambda)dlambda
	phi = probe_pulse_phase(wavelengths[j]);
	probe_pulse[i] = polar(A, phi);
    }
    for (; i<N_pulse; i++) probe_pulse[i] = 0;
    // calculate the peak intensity of the probe pulse
    pulse->Set(probe_pulse);
    pulse->ToTimeDomain();
    // make the peak intensity of the probe pulse equal to 1
    pulse->Normalise(Real(1));
    // prechirp the probe pulse
    pulse->ToFrequencyDomain();
    pulse->Get(probe_pulse);
    for (i=i1,j=0; j<number_of_frequencies; j++,i++)
    {
	phi = target_dispersion.Phase(frequencies[j]);
	probe_pulse[i] *= polar(Real(1), -phi);
    }
}

//-------------------------------------------------------------------------

PulseTarget::~PulseTarget()
{
   if (target_phase_shift) delete[] target_phase_shift;
   if (target_GDD) delete[] target_GDD;
   if (GDD) delete[] GDD;
   if (AD_weights) delete[] AD_weights;
   delete[] probe_pulse;
   delete[] reflected_pulse;
   delete pulse;
}

//-------------------------------------------------------------------------

Real PulseTarget::Merit(Coating& coating, Complex* reflectivity,
       Real* _reflectance, Real* _phase_shift)
{
    int i, j;
    Real I_max;
    int i1 = (N_pulse-number_of_frequencies)/2;

    if (!_reflectance) _reflectance = new Real[number_of_frequencies];
    if (!_phase_shift) _phase_shift = new Real[number_of_frequencies];
    
    if (!reflectivity) coating.Reflectivity(reflectivity);

    Reflectance(number_of_frequencies, reflectivity, _reflectance);
    PhaseShift(number_of_frequencies, reflectivity, _phase_shift);
    
    for (i=0; i<i1; i++) reflected_pulse[i] = 0;
    for (j=0; j<number_of_frequencies; j++)
	reflectivity[j] = polar(sqrt(_reflectance[j]), _phase_shift[j]);
    if (number_of_bounces > 1)
    {
	for (j=0; j<number_of_frequencies; j++)
	{
	    reflectivity[j] = IntPower(reflectivity[j], number_of_bounces);
	}
    }
    if (! adaptive_dispersion)
    {
	for (i=i1,j=0; j<number_of_frequencies; j++,i++)
	{
	    reflected_pulse[i] = probe_pulse[i] * reflectivity[j];
	}
    }
    else
    {
	int N = number_of_frequencies;
	Real omega_step = (frequencies[N-1] - frequencies[0]) / Real(N-1);
	PhaseToGDD(N, omega_step, _phase_shift, GDD);
	Real factor = AD_Factor(GDD);
	for (j=0; j<number_of_frequencies; j++,i++)
	{
	    reflected_pulse[i] = probe_pulse[i] * reflectivity[j] *
		polar(Real(1), target_phase_shift[j]*(1-factor));
	}
    } // if (! adaptive_dispersion)
    for (; i<N_pulse; i++) reflected_pulse[i] = 0;
    pulse->DestructiveDomainSwitch(false); // switch to the frequency domain
    pulse->Set(reflected_pulse);
    pulse->ToTimeDomain();
    I_max = pulse->PeakIntensity();
    return 1.0-I_max;
}

//-------------------------------------------------------------------------

Real PulseTarget::AD_Factor(const Real* GDD) const
{
    int i;
    Real s1=0, s2=0;
    if (! adaptive_dispersion) return 1;
    for (i=0; i<number_of_frequencies; i++)
    {
	s1 += target_GDD[i] * GDD[i] * AD_weights[i];
	s2 += square(target_GDD[i]) * AD_weights[i];
    }
    if (s2==Real(0)) return 1;
    return s1/s2;
}

//-------------------------------------------------------------------------
