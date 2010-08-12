#include "analysis.hh"
#include "coating.hh"
#include "target.hh"
#include "random.hh"
#include "pulse.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>

#define REPORT_FILE_NAME "analysis.txt"

const Real errorbar_ratio = 0.6827;

//======================================================================

void TechnologicalPerturbation(int no_of_layers, Real* orig_thicknesses,
			       Real* new_thicknesses,
			       Real mean_technological_error);

void SaveData(ostream& fs, int N, Real* X,
	      const vector< vector<Real> >& Y);
/* calculates the statistical properties (if necessary) and saves
   the data to the given file stream (which must be open);
   "X" is usually an array of wavelengths, "N" is its length, 
   the first dimension of "Y" is the perturbation index
   (zero means the unperturbed design), the second dimension of "Y"
   has the same meaning as "X" */

Real StatisticalAnalysis(vector<Real>::const_iterator start,
			 vector<Real>::const_iterator stop,
			 Real* standard_deviation=NULL);
/* returns the average value; if "standard_deviation" is not NULL,
   the standard deviation will be saved in this variable */


int ChoosePulseLength(int N, Real omega_step,
		      Real suggested_time_resolution);
// returns the number of points, which should be used for the probe pulse

void InitialiseProbePulse(int N_pulse, int N, Real omega_step,
			  AbstractParameter& probe_pulse_spectrum,
			  AbstractParameter& probe_pulse_phase,
			  AbstractDispersion& target_dispersion,
			  Real* wavelengths,
			  Pulse& pulse, Complex* probe_pulse);
/* Initialises the probe pulse. The array 'probe_pulse' must contain at least
   'N_pulse' elements, the array 'wavelengths' must contain at least 'N'
   elements */

void SaveProbeSpectrum(ostream& fs, int dots_in_plot,
		       const Real* wavelengths,
		       AbstractParameter& probe_pulse_spectrum,
		       AbstractParameter& probe_pulse_phase);
// saves the probe pulse in the frequency domain

void SaveProbePulse(ostream& fs, Pulse& pulse);
// saves the probe pulse in the time domain

//======================================================================

void Analyse(vector<Design>& designs, bool complementary_coatings,
	     const int dots_in_plot,
	     MaterialRepository& material_repository,
	     Parameters& parameters, Variables& variables)
{
    int i, j, k, no_of_layers, N_pulse=0;
    const int no_of_coatings = designs.size();
    int no_of_bounces, no_of_trials, merit_power;
    Real lambda_min, lambda_max, omega_min, omega_max;
    Real lambda0, omega0;
    Real mean_technological_error, individualism;
    Real omega_step, average, std, GD_contribution;
    Real target_EF, EF_tolerance, EF_reserve, dx;
    bool sensitivity_required, adaptive_dispersion, pulse_simulation;
    bool verbose, save_transmittance, save_absorption, lossy_materials;
    bool save_phase;
    Real wavelengths[dots_in_plot];
    Complex reflectivity[dots_in_plot];
    Real phase_shift[dots_in_plot];
    Real GD[dots_in_plot];
    Real target_GD[dots_in_plot];
    Real adapted_target_GD[dots_in_plot];
    Real target_GDD[dots_in_plot];
    Real target_phase_shift[dots_in_plot];
    Real GDD[dots_in_plot];
    Real GD_reserve[dots_in_plot];
    Real GD_tolerance[dots_in_plot];
    Real GDD_reserve[dots_in_plot];
    Real GDD_tolerance[dots_in_plot];
    Real reflectance[dots_in_plot];
    Real transmittance[dots_in_plot];
    Complex *probe_pulse = NULL;
    Complex *reflected_pulse = NULL;
    Real *pulse_intensity = NULL;
    Real *time = NULL;
    Pulse *pulse = NULL;
    vector< vector< vector<Complex> > > reflectivities; // for each design,
						       // each perturbation
						       // each wavelength
    vector< vector< vector<Real> > > property; // for each design,
                                               // each perturbation
                                               // each wavelength (time)
    vector< vector< vector<Real> > > transmittances; // for each design,
						     // each perturbation
						     // each wavelength
    vector< vector< Real > > merits; // for each design, for each perturbation
    vector<Real> adaptation_factor; // for each design
    string file_name, str;
    ofstream file;

    if (no_of_coatings==1) complementary_coatings = false;

    // get some parameters
    verbose << parameters("VERBOSE");
    no_of_bounces << parameters("NUMBER OF BOUNCES");
    no_of_trials << parameters("NUMBER OF TRIALS");
    lambda_min << parameters("MINIMAL WAVELENGTH");
    lambda_max << parameters("MAXIMAL WAVELENGTH");
    lambda0 << parameters("CENTRAL WAVELENGTH");
    omega0 = 2.0*M_PI*SPEED_OF_LIGHT / lambda0;
    omega_min = 2.0*M_PI*SPEED_OF_LIGHT / lambda_max;
    omega_max = 2.0*M_PI*SPEED_OF_LIGHT / lambda_min;
    mean_technological_error << parameters("TECHNOLOGICAL ERROR");
    sensitivity_required << parameters("SENSITIVITY");
    adaptive_dispersion << parameters("ADAPTIVE DISPERSION");
    save_phase << parameters("SAVE PHASE");
    if (complementary_coatings) individualism << parameters("INDIVIDUALISM");
    str << parameters("PROBE PULSE");
    pulse_simulation = (str.size() > 0);

    // create a vector of coatings
    Coating coating;
    vector<Coating> coatings;
    Coating::Parameters coating_parameters;
    coating_parameters.polarisation << parameters("POLARISATION");
    coating_parameters.angle_of_incidence << parameters("ANGLE OF INCIDENCE");
    coating_parameters.angle_of_incidence *= M_PI/180.0;
    coating_parameters.incidence_medium_name << parameters("INCIDENCE MEDIUM");
    coating_parameters.exit_medium_name << parameters("EXIT MEDIUM");
    coating_parameters.N = dots_in_plot;
    coating_parameters.omega_min = omega_min;
    coating_parameters.omega_max = omega_max;
    coating.SetParameters(coating_parameters, material_repository);
    for (i=0; i<no_of_coatings; i++) coatings.push_back(coating);
    for (i=0; i<no_of_coatings; i++)
	coatings[i].ImportDesign(designs[i], material_repository);

    // decide wether the transmittance and/or absorption have to be saved
    // (if not specified in the parameter file, find out whether some of the
    // materials are lossy (or have gain) in the wavelength range of interest
    lossy_materials = material_repository.ComplexIndex(lambda_min, lambda_max);
    str << parameters("SAVE TRANSMITTANCE");
    if (str=="") save_transmittance = lossy_materials;
    else save_transmittance << parameters("SAVE TRANSMITTANCE");
    str << parameters("SAVE ABSORPTION");
    if (str=="") save_absorption = lossy_materials;
    else save_absorption << parameters("SAVE ABSORPTION");
    lossy_materials = (save_transmittance || save_absorption);

    coating.Wavelengths(wavelengths);
    omega_step = (omega_max - omega_min) / Real(dots_in_plot-1);

    // reserve the necessary memory for the vectors
    if (! sensitivity_required) no_of_trials = 0;
    reflectivities.resize(no_of_coatings);
    property.resize(no_of_coatings);
    adaptation_factor.resize(no_of_coatings+1);
    for (i=0; i<no_of_coatings; i++)
    {
	reflectivities[i].resize(1+no_of_trials);
	property[i].resize(1+no_of_trials);
	for (j=0; j<=no_of_trials; j++)
	{
	    reflectivities[i][j].resize(dots_in_plot);
	    property[i][j].resize(dots_in_plot);
	}
    }
    if (lossy_materials)
    {
	transmittances.resize(no_of_coatings);
	for (i=0; i<no_of_coatings; i++)
	{
	    transmittances[i].resize(1+no_of_trials);
	    for (j=0; j<=no_of_trials; j++)
	    {
		transmittances[i][j].resize(dots_in_plot);
	    }
	}
    }
    if (complementary_coatings)
    {
	merits.resize(1);
	merits[0].resize(1+no_of_trials);
    }
    else
    {
	merits.resize(no_of_coatings);
	for (i=0; i<no_of_coatings; i++) merits[i].resize(1+no_of_trials);
    }

    for (i=0; i<=no_of_coatings; i++) adaptation_factor[i] = 1;

    // if the pulse analysis is requested, initialise the pulse
    // and allocate the memory
    if (pulse_simulation)
    {
	Real pulse_dt, t_min, time_resolution;
	time_resolution = 0.1 / (omega_max - omega_min);
	N_pulse = ChoosePulseLength(dots_in_plot, omega_step, 
				    time_resolution);
	pulse_dt = 2.0*M_PI / (N_pulse*omega_step);
	probe_pulse = new Complex[N_pulse];
	reflected_pulse = new Complex[N_pulse];
	pulse_intensity = new Real[N_pulse];
	time = new Real[N_pulse];
	t_min = -pulse_dt*(N_pulse/2);
	pulse = new Pulse(N_pulse, pulse_dt, t_min, false);
	for (i=0; i<N_pulse; i++) time[i] = t_min + i*pulse_dt;
	InitialiseProbePulse(N_pulse, dots_in_plot, omega_step,
			     variables.ProbePulseSpectrum(),
			     variables.ProbePulseSpectralPhase(),
			     variables.TargetDispersion(),
			     wavelengths, *pulse, probe_pulse);
	for (i=0; i<no_of_coatings; i++)
	{
	    for (j=0; j<=no_of_trials; j++) property[i][j].resize(N_pulse);
	}
    }

    // create a target
    AbstractTarget* target;
    string target_type;
    target_type << parameters("TARGET TYPE");
    if (target_type == "GDD")
    {
	GD_contribution << parameters("GD CONTRIBUTION");
	merit_power << parameters("MERIT POWER");
	
	target_EF << parameters("TARGET EF");
	EF_tolerance << parameters("EF TOLERANCE");
	EF_reserve << parameters("EF RESERVE");
	dx << parameters("STEP");

	target = new GDDTarget(dots_in_plot,
			       omega_min, omega_max,
			       variables.TargetDispersion(),
			       variables.TargetReflectance(),
			       variables.GDTolerance(),
			       variables.GDReserve(),
			       variables.GDDTolerance(),
			       variables.GDDReserve(),
			       GD_contribution,
			       variables.ReflectanceTolerance(),
			       variables.ReflectanceReserve(),
			       no_of_bounces,
			       target_EF,
			       EF_tolerance,
			       EF_reserve,
			       dx,
			       Real(-1), // no thickness constraint
			       individualism,
			       merit_power,
			       adaptive_dispersion);
    }
    else if (target_type == "pulse")
    {
	target = new PulseTarget(dots_in_plot,
				 omega_min, omega_max,
				 variables.ProbePulseSpectrum(),
				 variables.ProbePulseSpectralPhase(),
				 variables.TargetDispersion(),
				 0.2, // time resolution
				 no_of_bounces,
				 Real(-1), // no thickness constraint
				 individualism,
				 adaptive_dispersion);
    }
    else
    {
	cerr << "ERROR in 'Analyse': \
unknown optimisation target: " << target_type << "\n"
	     << "correct options: GDD, pulse." << endl;
	exit(1);
    }
    
    // calculate the target dispersion as well as the dispersion reserve
    // and tolerance
    {
	Real omega;
	AbstractDispersion& target_dispersion = variables.TargetDispersion();
	VariableParameter& GD_tolerance_parameter = variables.GDTolerance();
	VariableParameter& GD_reserve_parameter = variables.GDReserve();
	VariableParameter& GDD_tolerance_parameter = variables.GDDTolerance();
	VariableParameter& GDD_reserve_parameter = variables.GDDReserve();
	for (j=0; j<dots_in_plot; j++)
	{
	    omega = omega_min + j*omega_step;
	    target_GD[j] = target_dispersion.GD(omega);
	    target_GDD[j] = target_dispersion.GDD(omega);
	    target_phase_shift[j] = target_dispersion.Phase(omega);
	    GD_reserve[j] = GD_reserve_parameter(wavelengths[j]);
	    GD_tolerance[j] = GD_tolerance_parameter(wavelengths[j]);
	    GDD_reserve[j] = GDD_reserve_parameter(wavelengths[j]);
	    GDD_tolerance[j] = GDD_tolerance_parameter(wavelengths[j]);
    }
    }

    // analyse the designs
    for (i=0; i<no_of_coatings; i++)
    {
	coatings[i].Reflectivity(reflectivity);
	for (j=0; j<dots_in_plot; j++)
	    reflectivities[i][0][j] = reflectivity[j];
	if (! complementary_coatings)
	    merits[i][0] = target->MeritFunction(coatings[i]);
	if (lossy_materials)
	{
	    coatings[i].Transmittance(transmittance);
	    for (j=0; j<dots_in_plot; j++)
		transmittances[i][0][j] = transmittance[j];
	}
    }
    if (complementary_coatings) merits[0][0] = target->MeritFunction(coatings);
    // if the sensitivity analysis is required, simulate perturbations of
    // layer thicknesses
    if (sensitivity_required)
    {
	if (complementary_coatings)
	{
	    Real **thicknesses = new Real*[no_of_coatings];
	    Real **new_thicknesses = new Real*[no_of_coatings];
	    for (i=0; i<no_of_coatings; i++)
	    {
		no_of_layers = coatings[i].NumberOfLayers();
		thicknesses[i] = new Real[no_of_layers];
		new_thicknesses[i] = new Real[no_of_layers];
		coatings[i].GetLayerThicknesses(thicknesses[i]);
	    }
	    for (k=1; k<=no_of_trials; k++)
	    {
		for (i=0; i<no_of_coatings; i++)
		{
		    // perturb the thicknesses
		    no_of_layers = coatings[i].NumberOfLayers();
		    TechnologicalPerturbation(no_of_layers, thicknesses[i],
					      new_thicknesses[i],
					      mean_technological_error);
		    coatings[i].SetLayerThicknesses(new_thicknesses[i]);
		    // analyse the perturbed design
		    coatings[i].Reflectivity(reflectivity);
		    for (j=0; j<dots_in_plot; j++)
			reflectivities[i][k][j] = reflectivity[j];
		    if (lossy_materials)
		    {
			coatings[i].Transmittance(transmittance);
			for (j=0; j<dots_in_plot; j++)
			    transmittances[i][k][j] = transmittance[j];
		    }
		}
		merits[0][k] = target->MeritFunction(coatings);
	    }
	    for (i=0; i<no_of_coatings; i++)
	    {
		delete[] thicknesses[i];
		delete[] new_thicknesses[i];
	    }
	    delete[] thicknesses;
	    delete[] new_thicknesses;
	}
	else
	{
	    Real *thicknesses, *new_thicknesses;
	    for (i=0; i<no_of_coatings; i++)
	    {
		// perturb the thicknesses
		no_of_layers = coatings[i].NumberOfLayers();
		thicknesses = new Real[no_of_layers];
		new_thicknesses = new Real[no_of_layers];
		coatings[i].GetLayerThicknesses(thicknesses);
		for (k=1; k<=no_of_trials; k++)
		{
		    TechnologicalPerturbation(no_of_layers, thicknesses,
					      new_thicknesses,
					      mean_technological_error);
		    coatings[i].SetLayerThicknesses(new_thicknesses);
		    // analyse the perturbed design
		    coatings[i].Reflectivity(reflectivity);
		    for (j=0; j<dots_in_plot; j++)
			reflectivities[i][k][j] = reflectivity[j];
		    merits[i][k] = target->MeritFunction(coatings[i]);
		    if (lossy_materials)
		    {
			coatings[i].Transmittance(transmittance);
			for (j=0; j<dots_in_plot; j++)
			    transmittances[i][k][j] = transmittance[j];
		    }
		}
		delete[] thicknesses;
		delete[] new_thicknesses;
	    }
	} // if (complementary_coatings)
    } // if (sensitivity_required)

    // destroy the coatings
    coatings.clear();

    // analyse the reflectance
    for (i=0; i<no_of_coatings; i++)
    {
	for (k=0; k<=no_of_trials; k++)
	{
	    for (j=0; j<dots_in_plot; j++)
		reflectivity[j] = reflectivities[i][k][j];
	    Reflectance(dots_in_plot, reflectivity, reflectance);
	    if (no_of_bounces>1)
	    {
		for (j=0; j<dots_in_plot; j++)
		    reflectance[j] = IntPower(reflectance[j], no_of_bounces);
	    }
	    for (j=0; j<dots_in_plot; j++)
		property[i][k][j] = reflectance[j] * 100.0; // %
	}
	file_name = designs[i].name + ".R.dat";
	file.clear();
	file.open(file_name.c_str());
	file << "# The Reflectance of " << designs[i].name;
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: reflectance [%];\n";
	if (sensitivity_required)
	{
	    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	}
	SaveData(file, dots_in_plot, wavelengths, property[i]);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    } // for (i=0; i<no_of_coatings; i++)
    if (complementary_coatings)
    {   // calculate the averaged reflectance and its perturbations;
	// store them in property[0]
	for (k=0; k<=no_of_trials; k++)
	{
	    for (j=0; j<dots_in_plot; j++)
	    {
		for (i=0; i<no_of_coatings; i++)
		    property[i][k][j] /= 100.0;
		for (i=1; i<no_of_coatings; i++)
		    property[0][k][j] *= property[i][k][j];
		property[0][k][j] = 100.0 * pow(property[0][k][j],
						Real(1)/no_of_coatings);
	    }
	}
	file.clear();
	file_name = "averaged.R.dat";
	file.open(file_name.c_str());
	file << "# The averaged reflectance";
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: reflectance [%];\n";
	if (sensitivity_required)
	{
	    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	}
	SaveData(file, dots_in_plot, wavelengths, property[0]);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    } // if (complementary_coatings)


    // analyse the transmittance and the absorption
    if (lossy_materials)
    {
	for (i=0; i<no_of_coatings; i++)
	{
	    for (k=0; k<=no_of_trials; k++)
	    {
		if (no_of_bounces==1)
		{
		    for (j=0; j<dots_in_plot; j++)
			property[i][k][j] = transmittances[i][k][j] * 100;
		}
		else
		{
		    for (j=0; j<dots_in_plot; j++)
			reflectivity[j] = reflectivities[i][k][j];
		    Reflectance(dots_in_plot, reflectivity, reflectance);
		    for (j=0; j<dots_in_plot; j++)
		    {
			property[i][k][j] = 100 * transmittances[i][k][j] *
			    (1 - IntPower(reflectance[j], no_of_bounces)) /
			    (1 - reflectance[j]); // %
		    }
		}
	    } // for (k=0; k<=no_of_trials; k++)
	    if (save_transmittance)
	    {
		file_name = designs[i].name + ".T.dat";
		file.clear();
		file.open(file_name.c_str());
		file << "# The energy loss through Transmission in "
		     << designs[i].name;
		if (no_of_bounces>1) file << " (" << no_of_bounces
					  << " bounces)";
		file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: transmittance [%];\n";
		if (sensitivity_required)
		{
		    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
		}
		SaveData(file, dots_in_plot, wavelengths, property[i]);
		file.close();
		if (! file) cout << "ERROR: can't write " << file_name << endl;
		else if (verbose) cout << file_name << " has been saved"
				       << endl;
	    } // if (save_transmittance)
	    if (save_absorption)
	    {
		// calculate and save the absorption
		for (k=0; k<=no_of_trials; k++)
		{
		    for (j=0; j<dots_in_plot; j++)
			reflectivity[j] = reflectivities[i][k][j];
		    Reflectance(dots_in_plot, reflectivity, reflectance);
		    if (no_of_bounces>1)
		    {
			for (j=0; j<dots_in_plot; j++)
			{
			    reflectance[j] =
				IntPower(reflectance[j], no_of_bounces);
			}
		    }
		    for (j=0; j<dots_in_plot; j++)
			property[i][k][j] = 100 * (1 - reflectance[j]) -
			    property[i][k][j]; // %
		} // for (k=0; k<=no_of_trials; k++)
		file_name = designs[i].name + ".A.dat";
		file.clear();
		file.open(file_name.c_str());
		file << "# The energy loss through Absorption in "
		     << designs[i].name;
		if (no_of_bounces>1) file << " (" << no_of_bounces
					  << " bounces)";
		file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: absorption [%];\n";
		if (sensitivity_required)
		{
		    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
		}
		SaveData(file, dots_in_plot, wavelengths, property[i]);
		file.close();
		if (! file) cout << "ERROR: can't write " << file_name << endl;
		else if (verbose) cout << file_name << " has been saved"
				       << endl;
	    } // if (save_absorption)
	} // for (i=0; i<no_of_coatings; i++)
    } // if (lossy_materials)

    // analyse the GDD and calculate the adaptation factors
    for (i=0; i<no_of_coatings; i++)
    {
	for (k=0; k<=no_of_trials; k++)
	{
	    for (j=0; j<dots_in_plot; j++)
		reflectivity[j] = reflectivities[i][k][j];
	    PhaseShift(dots_in_plot, reflectivity, phase_shift);
	    PhaseToGDD(dots_in_plot, omega_step, phase_shift, GDD);
	    if (no_of_bounces>1)
		for (j=0; j<dots_in_plot; j++) GDD[j] *= no_of_bounces;
	    if (adaptive_dispersion && k==0)
	    {
		adaptation_factor[i] =
		    AdaptiveDispersionFactor(dots_in_plot, target_GDD,
					     GDD_reserve, GDD_tolerance, GDD);

	    }
	    for (j=0; j<dots_in_plot; j++) property[i][k][j] = GDD[j];
	}
	file_name = designs[i].name + ".GDD.dat";
	file.clear();
	file.open(file_name.c_str());
	file << "# The Group Delay Dispersion of " << designs[i].name;
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: GDD [fs^2];\n";
	if (sensitivity_required)
	{
	    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	}
	SaveData(file, dots_in_plot, wavelengths, property[i]);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    } // for (i=0; i<no_of_coatings; i++)
    if (complementary_coatings)
    {   // calculate the averaged GDD and its perturbations;
	// store them in property[0]
	for (k=0; k<=no_of_trials; k++)
	{
	    for (j=0; j<dots_in_plot; j++)
	    {
		for (i=1; i<no_of_coatings; i++)
		    property[0][k][j] += property[i][k][j];
		property[0][k][j] /= Real(no_of_coatings);
	    }
	    if (adaptive_dispersion && k==0)
	    {
		for (j=0; j<dots_in_plot; j++) GDD[j] = property[0][k][j];
		adaptation_factor[no_of_coatings] =
		    AdaptiveDispersionFactor(dots_in_plot, target_GDD,
					     GDD_reserve, GDD_tolerance, GDD);
	    }
	}
	file.clear();
	file_name = "averaged.GDD.dat";
	file.open(file_name.c_str());
	file << "# The averaged Group Delay Dispersion";
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: GDD [fs^2];\n";
	if (sensitivity_required)
	{
	    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	}
	SaveData(file, dots_in_plot, wavelengths, property[0]);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    } // if (complementary_coatings)

    // analyse the phase shift upon reflection
    if (save_phase)
    {
	for (i=0; i<no_of_coatings; i++)
	{
	    for (k=0; k<=no_of_trials; k++)
	    {
		for (j=0; j<dots_in_plot; j++)
		    reflectivity[j] = reflectivities[i][k][j];
		PhaseShift(dots_in_plot, reflectivity, phase_shift);
		if (no_of_bounces>1)
		{
		    for (j=0; j<dots_in_plot; j++)
			phase_shift[j] *= no_of_bounces;
		}
		for (j=0; j<dots_in_plot; j++)
		    property[i][k][j] = phase_shift[j];
	    }
	    file_name = designs[i].name + ".phase.dat";
	    file.clear();
	    file.open(file_name.c_str());
	    file << "# The phase shift upon reflection off "
		 << designs[i].name;
	    if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	    file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: phase shift [radians];\n";
	    if (sensitivity_required)
	    {
		file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	    }
	    SaveData(file, dots_in_plot, wavelengths, property[i]);
	    file.close();
	    if (! file) cout << "ERROR: can't write " << file_name << endl;
	    else if (verbose) cout << file_name << " has been saved" << endl;
	} // for (i=0; i<no_of_coatings; i++)
	if (complementary_coatings)
	{   // calculate the averaged phase shift and its perturbations;
	    // store them in property[0]
	    for (k=0; k<=no_of_trials; k++)
	    {
		for (j=0; j<dots_in_plot; j++)
		{
		    for (i=1; i<no_of_coatings; i++)
			property[0][k][j] += property[i][k][j];
		    property[0][k][j] /= Real(no_of_coatings);
		}
	    }
	    file.clear();
	    file_name = "averaged.phase.dat";
	    file.open(file_name.c_str());
	    file << "# The averaged phase shift upon reflection";
	    if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	    file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: phase shift [radians];\n";
	    if (sensitivity_required)
	    {
		file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	    }
	    SaveData(file, dots_in_plot, wavelengths, property[0]);
	    file.close();
	    if (! file) cout << "ERROR: can't write " << file_name << endl;
	    else if (verbose) cout << file_name << " has been saved" << endl;
	} // if (complementary_coatings)
    } // if (save_phase)
    
    // analyse the GD
    for (i=0; i<no_of_coatings; i++)
    {
	for (k=0; k<=no_of_trials; k++)
	{
	    for (j=0; j<dots_in_plot; j++)
		reflectivity[j] = reflectivities[i][k][j];
	    PhaseShift(dots_in_plot, reflectivity, phase_shift);
	    PhaseToGD(dots_in_plot, omega_step, phase_shift, GD);
	    if (no_of_bounces>1)
		for (j=0; j<dots_in_plot; j++) GD[j] *= no_of_bounces;
	    for (j=0; j<dots_in_plot; j++) property[i][k][j] = GD[j];
	}
	file_name = designs[i].name + ".GD.dat";
	file.clear();
	file.open(file_name.c_str());
	file << "# The Group Delay of " << designs[i].name;
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: GD [fs];\n";
	if (sensitivity_required)
	{
	    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	}
	SaveData(file, dots_in_plot, wavelengths, property[i]);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    } // for (i=0; i<no_of_coatings; i++)
    if (complementary_coatings)
    {   // calculate the averaged GD and its perturbations;
	// store them in property[0]
	for (k=0; k<=no_of_trials; k++)
	{
	    for (j=0; j<dots_in_plot; j++)
	    {
		for (i=1; i<no_of_coatings; i++)
		    property[0][k][j] += property[i][k][j];
		property[0][k][j] /= Real(no_of_coatings);
	    }
	}
	file.clear();
	file_name = "averaged.GD.dat";
	file.open(file_name.c_str());
	file << "# The averaged Group Delay";
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: GD [fs];\n";
	if (sensitivity_required)
	{
	    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	}
	SaveData(file, dots_in_plot, wavelengths, property[0]);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    } // if (complementary_coatings)

    // analyse the uncompensated GD
    for (i=0; i<no_of_coatings; i++)
    {
	Real GD_shift;
	for (j=0; j<dots_in_plot; j++)
	    adapted_target_GD[j] = target_GD[j] * adaptation_factor[i];
	for (k=0; k<=no_of_trials; k++)
	{
	    for (j=0; j<dots_in_plot; j++)
		reflectivity[j] = reflectivities[i][k][j];
	    PhaseShift(dots_in_plot, reflectivity, phase_shift);
	    PhaseToGD(dots_in_plot, omega_step, phase_shift, GD);
	    if (no_of_bounces>1)
		for (j=0; j<dots_in_plot; j++) GD[j] *= no_of_bounces;
	    GD_shift = GDShift(dots_in_plot, adapted_target_GD, GD_reserve,
			       GD_tolerance, GD);
	    for (j=0; j<dots_in_plot; j++)
	    {
		property[i][k][j] = GD[j] + GD_shift - adapted_target_GD[j];
	    }
	}
	file_name = designs[i].name + ".uGD.dat";
	file.clear();
	file.open(file_name.c_str());
	file << "# The uncompensated Group Delay of " << designs[i].name;
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: GD [fs];\n";
	if (sensitivity_required)
	{
	    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	}
	SaveData(file, dots_in_plot, wavelengths, property[i]);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    } // for (i=0; i<no_of_coatings; i++)
    if (complementary_coatings)
    {   // calculate the averaged uncompensated GD and its perturbations;
	// store them in property[0]
	for (k=0; k<=no_of_trials; k++)
	{
	    for (j=0; j<dots_in_plot; j++)
	    {
		for (i=1; i<no_of_coatings; i++)
		    property[0][k][j] += property[i][k][j];
		property[0][k][j] /= Real(no_of_coatings);
	    }
	}
	file.clear();
	file_name = "averaged.uGD.dat";
	file.open(file_name.c_str());
	file << "# The averaged uncompensated Group Delay";
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: wavelength [nm];\n\
# column 2: GD [fs];\n";
	if (sensitivity_required)
	{
	    file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	}
	SaveData(file, dots_in_plot, wavelengths, property[0]);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    } // if (complementary_coatings)

    // if necessary, analyse the pulse reflection
    if (pulse_simulation)
    {
	int m;
	int m1 = (N_pulse-dots_in_plot)/2;
	file_name = "probe_spectrum.dat";
	file.clear();
	file.open(file_name.c_str());
	SaveProbeSpectrum(file, dots_in_plot, wavelengths,
			  variables.ProbePulseSpectrum(),
			  variables.ProbePulseSpectralPhase());
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
	for (m=0; m<m1; m++) reflected_pulse[m] = 0;
	for (m=m1,j=0; j<dots_in_plot; j++,m++)
	    reflected_pulse[m] = probe_pulse[m];
	for (; m<N_pulse; m++) reflected_pulse[m] = 0;
	pulse->DestructiveDomainSwitch(false);
	pulse->Set(reflected_pulse);
	pulse->ToTimeDomain();
	pulse->ShiftToCentre(omega0);
	pulse->Get(reflected_pulse);
	file_name = "probe_pulse.dat";
	file.clear();
	file.open(file_name.c_str());
	SaveProbePulse(file, *pulse);
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
	for (i=0; i<no_of_coatings; i++)
	{
	    for (k=0; k<=no_of_trials; k++)
	    {
		for (m=0; m<m1; m++) reflected_pulse[m] = 0;
		for (j=0; j<dots_in_plot; j++)
		    reflectivity[j] = reflectivities[i][k][j];
		if (no_of_bounces > 1)
		{
		    for (j=0; j<dots_in_plot; j++)
		    {
			reflectivity[j] =
			    IntPower(reflectivity[j], no_of_bounces);
		    }
		}
		for (m=m1,j=0; j<dots_in_plot; j++,m++)
		{
		    reflected_pulse[m] = probe_pulse[m] * reflectivity[j] *
			polar(Real(1),
			      -target_phase_shift[j]*adaptation_factor[i]);
		}
		for (; m<N_pulse; m++) reflected_pulse[m] = 0;
		pulse->DestructiveDomainSwitch(false);
		pulse->Set(reflected_pulse);
		pulse->ToTimeDomain();
		pulse->ShiftToCentre(omega0);
		pulse->GetIntensity(pulse_intensity);
		for (m=0; m<N_pulse; m++)
		    property[i][k][m] = pulse_intensity[m];
	    } // for (k=0; k<=no_of_trials; k++)
	    file_name = designs[i].name + ".I.dat";
	    file.clear();
	    file.open(file_name.c_str());
	    file << "\
# The probe pulse prechirped with the target dispersion and reflected off\n\
# the coating '" << designs[i].name << "'";
	    if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	    file << ":\n\
# column 1: time [fs];\n\
# column 2: intensity [arbitrary units];\n";
	    if (sensitivity_required)
	    {
		file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	    }
	    SaveData(file, N_pulse, time, property[i]);
	    file.close();
	    if (! file) cout << "ERROR: can't write " << file_name << endl;
	    else if (verbose) cout << file_name << " has been saved" << endl;
	} // for (i=0; i<no_of_coatings; i++)
	if (complementary_coatings)
	{   // calculate the pulse reflected off the set of mirrors
	    Real* net_reflectance = new Real[dots_in_plot];
	    Real* net_phase_shift = new Real[dots_in_plot];
	    for (k=0; k<=no_of_trials; k++)
	    {
		for (j=0; j<dots_in_plot; j++)
		{
		    net_reflectance[j] = 1;
		    net_phase_shift[j] = 0;
		}
		for (i=0; i<no_of_coatings; i++)
		{
		    for (j=0; j<dots_in_plot; j++)
			reflectivity[j] = reflectivities[i][k][j];
		    Reflectance(dots_in_plot, reflectivity, reflectance);
		    PhaseShift(dots_in_plot, reflectivity, phase_shift);
		    for (j=0; j<dots_in_plot; j++)
		    {
			net_reflectance[j] *= reflectance[j];
			net_phase_shift[j] += phase_shift[j];
		    }
		}
		for (j=0; j<dots_in_plot; j++)
		{
		    net_reflectance[j] =
			pow(net_reflectance[j], 1/Real(no_of_coatings));
		    net_phase_shift[j] /= no_of_coatings;
		}
		for (j=0; j<dots_in_plot; j++)
		{
		    reflectivity[j] = 1;
		    for (i=0; i<no_of_coatings; i++)
			reflectivity[j] *= reflectivities[i][k][j];
		}
		if (no_of_bounces > 1)
		{
		    for (j=0; j<dots_in_plot; j++)
		    {
			reflectivity[j] =
			    IntPower(reflectivity[j], no_of_bounces);
			net_reflectance[j] =
			    IntPower(net_reflectance[j], no_of_bounces);
			net_phase_shift[j] *= no_of_bounces;
		    }
		}
		// analyse the pulse reflected off the set of mirrors
		for (m=0; m<m1; m++) reflected_pulse[m] = 0;
		for (m=m1,j=0; j<dots_in_plot; j++,m++)
		{
		    reflected_pulse[m] = probe_pulse[m] * reflectivity[j] *
			polar(Real(1), -no_of_coatings*target_phase_shift[j] *
			      adaptation_factor[no_of_coatings]);
		}
		for (; m<N_pulse; m++) reflected_pulse[m] = 0;
		pulse->DestructiveDomainSwitch(false);
		pulse->Set(reflected_pulse);
		pulse->ToTimeDomain();
		pulse->ShiftToCentre(omega0);
		pulse->GetIntensity(pulse_intensity);
		for (m=0; m<N_pulse; m++)
		    property[0][k][m] = pulse_intensity[m];
		// analyse the pulse after one "virtual" bounce
		for (j=0; j<dots_in_plot; j++)
		{
		    reflectivity[j] = polar(sqrt(net_reflectance[j]),
						 net_phase_shift[j]);
		}
		for (m=0; m<m1; m++) reflected_pulse[m] = 0;
		for (m=m1,j=0; j<dots_in_plot; j++,m++)
		{
		    reflected_pulse[m] = probe_pulse[m] * reflectivity[j] *
			polar(Real(1), -target_phase_shift[j] *
			      adaptation_factor[no_of_coatings]);
		}
		for (; m<N_pulse; m++) reflected_pulse[m] = 0;
		pulse->DestructiveDomainSwitch(false);
		pulse->Set(reflected_pulse);
		pulse->ToTimeDomain();
		pulse->ShiftToCentre(omega0);
		pulse->GetIntensity(pulse_intensity);
		for (m=0; m<N_pulse; m++)
		    property[1][k][m] = pulse_intensity[m];
	    } // for (k=0; k<=no_of_trials; k++)
	    file.clear();
	    file_name = "altogether.I.dat";
	    file.open(file_name.c_str());
	    file << "\
# The probe pulse prechirped with the target dispersion and reflected off\n\
# the combination of the coatings";
	    if (no_of_bounces>1)
	    {
		file << " (" << no_of_bounces << " bounces off each of \
the coatings)";
	    }
	    file << ":\n\
# column 1: time [fs];\n\
# column 2: intensity [arbitrary units];\n";
	    if (sensitivity_required)
	    {
		file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	    }
	    SaveData(file, N_pulse, time, property[0]);
	    file.close();
	    if (! file) cout << "ERROR: can't write " << file_name << endl;
	    else if (verbose) cout << file_name << " has been saved" << endl;
	    file.clear();
	    file_name = "averaged.I.dat";
	    file.open(file_name.c_str());
	    file << "\
# The probe pulse prechirped with the target dispersion and reflected off\n\
# a \"virtual\" mirror, which represents the averaged properties of the whole\n\
# system of mirrors";
	if (no_of_bounces>1) file << " (" << no_of_bounces << " bounces)";
	file << ":\n\
# column 1: time [fs];\n\
# column 2: intensity [arbitrary units];\n";
	    if (sensitivity_required)
	    {
		file << "# column 3: the minimal value of the error bar;\n\
# column 4: the maximal value of the error bar.\n";
	    }
	    SaveData(file, N_pulse, time, property[1]);
	    file.close();
	    if (! file) cout << "ERROR: can't write " << file_name << endl;
	    else if (verbose) cout << file_name << " has been saved" << endl;
	    delete[] net_reflectance;
	    delete[] net_phase_shift;
	} // if (complementary_coatings)
    } // if (pulse_simulation)

    if (adaptive_dispersion)
    {   // save the adaptation factors
	file.clear();
	file_name = "adaptation.txt";
	file.open(file_name.c_str());
	file << "\
# The analysis was made with the ADAPTIVE DISPERSION turned on. This file\n\
# contains the factors, by which the target dispersion was multiplied.\n";
	for (i=0; i<no_of_coatings; i++)
	    file << designs[i].name << "\t" << adaptation_factor[i] << "\n";
	if (complementary_coatings)
	    file << "averaged\t" << adaptation_factor[no_of_coatings] << "\n";
	file.close();
	if (! file) cout << "ERROR: can't write " << file_name << endl;
	else if (verbose) cout << file_name << " has been saved" << endl;
    }

    // save the information about the merit function
    file.clear();
    file_name = REPORT_FILE_NAME;
    file.open(file_name.c_str());
    if (complementary_coatings)
	file << "# The merit function of the set of the mirrors\n";
    file << "# Here is the merit function of the analysed designs.\n";
    if ( sensitivity_required )
    {
	file << "\
# This file also contains some statistical information about what random\n\
# technological perturbations can do with the merit function.\n\
# Averaging was done over "<< no_of_trials <<" randomly perturbed designs.\n\
# The mean technological error was taken to be equal to "<< 
	    mean_technological_error <<" nm.\n";
    }
    file <<  "# column 1: design name;\n\
# column 2: the figure of merit;\n";
    if ( sensitivity_required )
    {
	file << "\
# column 3: the mean value of the merit function;\n\
# column 4: the standard deviation of the merit function.\n";
    }
    if (complementary_coatings)
    {
	file << "averaged " << merits[0][0];
	if (sensitivity_required)
	{
	    vector<Real>::iterator it = merits[0].begin();
	    it++;
	    average = StatisticalAnalysis(it, merits[0].end(), &std);
	    file << " " << setprecision(4) << average <<" " << std << "\n";
	}
	else file << "\n";
    }
    else
    {
	for (i=0; i<no_of_coatings; i++)
	{
	    file << designs[i].name << " " << merits[i][0];
	    if (sensitivity_required)
	    {
		vector<Real>::iterator it = merits[i].begin();
		it++;
		average = StatisticalAnalysis(it, merits[i].end(), &std);
		file << " " << setprecision(4) << average <<" " << std << "\n";
	    }
	    else file << "\n";
	}
    }
    file.close();
    if (! file) cout << "ERROR: can't write " << file_name << endl;
    else if (verbose) cout << file_name << " has been saved" << endl;
    delete target;
    if (probe_pulse) delete[] probe_pulse;
    if (reflected_pulse) delete[] reflected_pulse;
    if (pulse_intensity) delete[] pulse_intensity;
    if (time) delete[] time;
    if (pulse) delete pulse;
}

//=========================================================================

void TechnologicalPerturbation(int no_of_layers, Real* orig_thicknesses,
			       Real* new_thicknesses,
			       Real mean_technological_error)
{
    Real random_shift;
    int layer;
    for (layer=0; layer<no_of_layers; layer++)
    {
	if (orig_thicknesses[layer] < 0)
	    orig_thicknesses[layer] = - orig_thicknesses[layer];
	do
	{
	    random_shift = mean_technological_error * GaussianRandom();
	    new_thicknesses[layer] = orig_thicknesses[layer] + random_shift;
	} while (new_thicknesses[layer] <= 0);
    }
}

//-------------------------------------------------------------------------

void CalculateErrorbars(vector<Real>& data,
			Real *errorbar_min, Real *errorbar_max)
{
    int N = data.size();
    int j, n;
    int bar_length = (int) round((Real)N * errorbar_ratio);
    Real length, min_length;
    sort(data.begin(), data.end());
    // determine the position of the error bar
    n = 0;
    min_length = data[bar_length] - data[0];
    for (j=1; j<N-bar_length; j++)
    {
	length = data[j+bar_length] - data[j];
	if (min_length > length)
	{
	    min_length = length;
	    n = j;
	}
    }
    *errorbar_min = data[n];
    *errorbar_max = data[n+bar_length];
}

//-------------------------------------------------------------------------

void SaveData(ostream& fs, int N, Real* X,
	      const vector< vector<Real> >& Y)
{
    int i, j;
    int no_of_trials = Y.size() - 1;
    bool statistics = (no_of_trials > 0);
    Real errorbar_min, errorbar_max;
    vector<Real> data(no_of_trials);
    for (j=0; j<N; j++)
    {
	fs << fixed << setprecision(5) << setw(15) << X[j] << " "
	   << scientific << setprecision(7) << setw(17) << Y[0][j];
	if (! statistics) fs << "\n";
	else
	{
	    // calculate the position of the error bar
	    for (i=0; i<no_of_trials; i++) data[i] = Y[i+1][j];
	    CalculateErrorbars(data, &errorbar_min, &errorbar_max);
	    // save the position of the error bar
	    fs << " " << scientific << setprecision(7) << setw(15)
	       << errorbar_min << " " << errorbar_max << "\n";
	}
    }
}

//-------------------------------------------------------------------------

Real StatisticalAnalysis(vector<Real>::const_iterator start,
			 vector<Real>::const_iterator stop,
			 Real* standard_deviation)
{
    int count = 0;
    Real average = 0;
    vector<Real>::const_iterator it;
    for (it=start; it!=stop; it++)
    {
	average += *it;
	count++;
    }
    if (! count) return 0;
    average /= count;
    if (standard_deviation)
    {
	*standard_deviation = 0;
	for (it=start; it!=stop; it++)
	    *standard_deviation += square(*it - average);
	*standard_deviation /= count;
    }
    return average;
}

//-------------------------------------------------------------------------

int ChoosePulseLength(int N, Real omega_step,
		      Real suggested_time_resolution)
{
    int N_pulse;
    Real pulse_dt;
    // determine the number of points for the probe pulse
    N_pulse = IntPower(2, (int) ceil(log(Real(N))/ log(Real(2))));
    pulse_dt = 2.0*M_PI / (N_pulse*omega_step);
    if (suggested_time_resolution <= Real(0))
    {
	throw "ERROR in ChoosePulseLength: \
suggested_time_resolution <= 0";
    }
    while (pulse_dt > suggested_time_resolution*1.1)
    {
	N_pulse *= 2;
	pulse_dt /= 2.0;
    }
    return N_pulse;
}

//-------------------------------------------------------------------------

void InitialiseProbePulse(int N_pulse, int N, Real omega_step,
			  AbstractParameter& probe_pulse_spectrum,
			  AbstractParameter& probe_pulse_phase,
			  AbstractDispersion& target_dispersion,
			  Real* wavelengths,
			  Pulse& pulse, Complex* probe_pulse)
{
    int i, i1, j;
    Real A, phi;
    Real* frequencies = new Real[N];

    for (i=0; i<N; i++)
	frequencies[i] = 2.0*M_PI*SPEED_OF_LIGHT / wavelengths[i];

    // calculate the probe pulse (do not prechirp)
    if (N_pulse < N) throw "ERROR in InitialiseProbePulse: N_pulse < N";
    i1 = (N_pulse-N)/2;
    for (i=0; i<i1; i++) probe_pulse[i] = 0;
    for (i=i1,j=0; j<N; j++,i++)
    {
	A = sqrt(probe_pulse_spectrum(wavelengths[j]));
	A /= frequencies[j]; // I(omega)domega = I(lambda)dlambda
	phi = probe_pulse_phase(wavelengths[j]);
	probe_pulse[i] = polar(A, phi);
    }
    for (; i<N_pulse; i++) probe_pulse[i] = 0;
    // calculate the peak intensity of the probe pulse
    pulse.ToFrequencyDomain();
    pulse.Set(probe_pulse);
    pulse.ToTimeDomain();
    // make the peak intensity of the probe pulse equal to 1
    pulse.Normalise(Real(1));
    // save the probe pulse in 'probe_pulse'
    pulse.ToFrequencyDomain();
    pulse.Get(probe_pulse);
    delete[] frequencies;
}

//-------------------------------------------------------------------------

void SaveProbeSpectrum(ostream& fs, int N,
		       const Real* wavelengths,
		       AbstractParameter& probe_pulse_spectrum,
		       AbstractParameter& probe_pulse_phase)
{
    int i;
    Real w;
    fs << "\
# The probe pulse in the frequency domain:\n\
# column 1: wavelength [nm];\n\
# column 2: spectral intensity [arbitrary units];\n\
# column 3: spectral phase [radians].\n";
    for (i=0; i<N; i++)
    {
	w = wavelengths[i];
	fs << fixed << setprecision(5) << setw(15) << w << " "
	   << scientific << setprecision(7) << setw(17)
	   << probe_pulse_spectrum(w) << " "
	   << setprecision(7) << setw(17) << probe_pulse_phase(w) << "\n";
    }
}

//-------------------------------------------------------------------------

void SaveProbePulse(ostream& fs, Pulse& pulse)
{
    int i, i1, i2;
    Real threshold, t, t_min, dt;
    int N = pulse.NumberOfPoints();
    Real *intensity = new Real[N];
    Real *phase = new Real[N];
    fs << "\
# The probe pulse in the time domain:\n\
# column 1: time [fs];\n\
# column 2: intensity [arbitrary units];\n\
# column 3: phase [radians].\n";
    pulse.ToTimeDomain();
    pulse.GetIntensity(intensity);
    pulse.GetPhase(phase);
    // determine the region, which we are going to save
    threshold = 0.001 * pulse.PeakIntensity();
    for (i1=0; i1<N && intensity[i1]<=threshold; i1++) {}
    if (i1==N) i1 = 0;
    for (i2=N-1; i2>=0 && intensity[i2]<=threshold; i2--) {}
    if (i2==0) i2 = N-1;
    // save the pulse in this region
    t_min = pulse.StartingTime();
    dt = pulse.TimeStep();
    for (i=i1; i<=i2; i++)
    {
	t = t_min + i*dt;
	fs << fixed << setprecision(5) << setw(15) << t << " "
	   << scientific << setprecision(7) << setw(17) << intensity[i]
	   << " " << setprecision(7) << setw(17) << phase[i] << "\n";
    }
    delete[] intensity;
    delete[] phase;
}

//-------------------------------------------------------------------------
