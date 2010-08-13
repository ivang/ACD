// Construct an initial design assuming that layer thickness varies smoothly

#include "common.hh"
#include "parameters.hh"
#include "material.hh"
#include "design.hh"
#include "coating.hh"

#include "target.hh"

// #include "optimisation.hh"
// #include "analysis.hh"
// #include "random.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>

#include "gsl/gsl_multimin.h"

//-------------------------------------------------------------------------

void Usage(Parameters& parameters)
{
    string s;
    cout << "Construction of an initial chirped-mirror design.\n\n\
USAGE: acdinit [options] <template designs>\n" << endl;
    s = "The template designs specify the materials and the number of \
layers in the stack, while all layer thicknesses in these files are \
ignored. Constructed designs are saved under the same names as the \
template ones in the current working directory.";
    cout << "  ";
    PrintFormatted(cout, s, 75, "");
    parameters.PrintHelp();
}

void InitialiseParameters(Parameters& parameters);

void InitialiseUnusedParameters(Parameters& parameters);
/* just to suppress the warnings */

void DesignNames(int argc, const char** argv, vector<string>& names);
// extracts the design names from the command line

void ConstructDesigns(vector<Design>& designs, Parameters& parameters,
		      Variables& variables,
		      MaterialRepository& material_repository,
		      AbstractTarget& target);

//-------------------------------------------------------------------------

class AbstractParameterConverter
{ /* Classes inherited from this class can convert a vector of abstract
   parameters into layer thicknesses, which allows optimisation of
   coatings, in which the thicknesses are given by some functional
   dependence. */
public:
    AbstractParameterConverter();
    virtual ~AbstractParameterConverter();
    virtual void LayerThicknesses(int number_of_parameters,
				  Real* parameters,
				  int number_of_layers,
				  Real *layer_thicknesses) = 0;
    void UpdateCoating(int number_of_parameters, Real* parameters,
		       Coating& coating);
private:
    int n_layers;
    Real* thicknesses;
};

class ChebyshevParameterConverter : public AbstractParameterConverter
{
public:
    ChebyshevParameterConverter();
    void LayerThicknesses(int number_of_parameters,
			  Real* parameters,
			  int number_of_layers,
			  Real *layer_thicknesses);
};

//-------------------------------------------------------------------------

void ParametricCGoptimisation(int number_of_parameters,
			      Real* parameters,
			      AbstractParameterConverter& converter,
			      Coating& coating, AbstractTarget& target,
			      Real step);
// conjugate-gradient optimisation

void ParametricSimplexOptimisation(int number_of_parameters,
				   Real* parameters,
				   AbstractParameterConverter& converter,
				   Coating& coating, AbstractTarget& target,
				   Real step);
// simplex optimisation; the 'step' here has nothing to do with the
// 'SIMPLEX STEP' in the parameter file, because our optimisation parameters
// are not layer thicknesses

//-------------------------------------------------------------------------

// the following global variables are initialised when a local optimisation
// starts

AbstractTarget* global_target = NULL;
Coating* optimised_coating = NULL;
AbstractParameterConverter* global_parameter_converter = NULL;
Real* global_optimisation_parameters = NULL;
int global_n_parameters = 0;
float global_gradient_step = 0;

//-------------------------------------------------------------------------

int main(int argc, const char** argv)
{
    int i, number_of_designs;
    vector<string> design_names;
    string str, str1, target_type;
    AbstractTarget* target = NULL;
    int number_of_frequencies, number_of_bounces, merit_power;
    Real lambda0, lambda_min, lambda_max, omega0, omega_min, omega_max;
    Real minimal_thickness;
    bool adaptive_dispersion;
    Real GD_contribution;
    Real target_EF, EF_tolerance, EF_reserve, dx;

    // initalise the parameters
    Parameters parameters;
    InitialiseParameters(parameters);
    if (argc<2) Usage(parameters);
    str = argv[1];
    if (str=="-h" || str=="--help") Usage(parameters);

    // read the parameter file
    InitialiseUnusedParameters(parameters); // to suppress the warnings
    {
	ifstream file("parameters.txt");
	if (! file)
	{
	    ofstream os("parameters.txt");
	    cerr << "Can't open the file 'parameters.txt'.\n\
Writing a template of this file in the current directory." << endl;
	    parameters.WriteParameterFile(os);
	    os.close();
	    cerr << "Edit this file and run the program again." << endl;
	    return 1;
	}
	try
	{
	    parameters.ParseFile(file);
	}
	catch (const string& s)
	{
	    cerr << "FAILED to parse the parameter file:\n"
		 << s << endl;
	    return 1;
	}
	file.close();
    }
    // parse the command line
    try
    {
	parameters.ParseCommandLine(argc-1, &argv[1]);
    }
    catch (const string& s)
    {
	cerr << "FAILED to parse the command line:\n"
	     << s << endl;
	return 1;
    }

    // get some parameters
    try
    {
	number_of_frequencies << parameters("RESOLUTION");
	lambda0 << parameters("CENTRAL WAVELENGTH");
	lambda_min << parameters("MINIMAL WAVELENGTH");
	lambda_max << parameters("MAXIMAL WAVELENGTH");
	omega0 = 2.0*M_PI*SPEED_OF_LIGHT / lambda0;
	omega_min = 2.0*M_PI*SPEED_OF_LIGHT / lambda_max;
	omega_max = 2.0*M_PI*SPEED_OF_LIGHT / lambda_min;
	number_of_bounces << parameters("NUMBER OF BOUNCES");
	minimal_thickness << parameters("MINIMAL THICKNESS");
	adaptive_dispersion << parameters("ADAPTIVE DISPERSION");
    }
    catch (const string& s)
    {
	cerr << "FAILED to get one of the important parameters:\n"
	     << s << endl;
	return 1;
    }

    // there are a few parameters, which are allowed to be uninitialised
    // at this point; we must initialise them now
    str << parameters("PULSE CENTRAL WAVELENGTH");
    if (str.size() == 0)
    {
	str << parameters("CENTRAL WAVELENGTH");
	parameters.SetParameterValue("PULSE CENTRAL WAVELENGTH", str.c_str());
    }
    str << parameters("PROBE PULSE FWHM");
    if (str.size() == 0)
    {
	int m;
	Real x, FWHM;
	ostringstream os;
	str << parameters("PROBE PULSE");
	istringstream is(str);
	is >> str;
	transform(str.begin(), str.end(), str.begin(), (int(*)(int))toupper);
	if (str!="SUPERGAUSSIAN") m = 1;
	else
	{
	    is >> m;
	    if (! is) m = 1;
	}
	if (lambda0-lambda_min < lambda_max-lambda0)
	    x = lambda0 - lambda_min;
	else x = lambda_max - lambda0;
	FWHM = 2.0*x * pow(log(2.0)/log(100.0), 0.5/m);
	os << FWHM;
	str = os.str();
	parameters.SetParameterValue("PROBE PULSE FWHM", str.c_str());
    }

    Variables variables(parameters);

    // read the design names
    DesignNames(argc-1, &argv[1], design_names);

    if (design_names.size() == 0)
    {
	cerr << "\
You must give one or more names of design files on the command line."
	     << endl;
	return 1;
    }

    try
    {
	// create the material repository
	vector<string> paths_to_materials;
	str << parameters("PATH TO MATERIALS");
	paths_to_materials.push_back(str);
#ifdef MS_WINDOWS
	// substitute '/' with '\' in 'PATH TO MATERIALS'
	str1 = str;
	replace(str1.begin(), str1.end(), '/', '\\');
	if (str != str1) paths_to_materials.push_back(str1);
#else
	// substitute '\' with '/' in 'PATH TO MATERIALS'
	str1 = str;
	replace(str1.begin(), str1.end(), '\\', '/');
	if (str != str1) paths_to_materials.push_back(str1);
#endif
	str = DEFAULT_PATH_TO_MATERIALS;
	paths_to_materials.push_back(str);
	MaterialRepository material_repository(paths_to_materials);
#ifdef DEBUG
	cout << "paths_to_materials:\n";
	for (i=0; i<paths_to_materials.size(); i++)
	    cout << "\t" << paths_to_materials[i] << endl;
	exit(0);
#endif

	// create a vector of designs
	Design design;
	vector<Design> designs;
	number_of_designs = design_names.size();
	for (i=0; i<number_of_designs; i++)
	{
	    design.ReadFromFile(design_names[i]);
	    designs.push_back(design);
	}

	// create a target
	target_type << parameters("TARGET TYPE");
	if (target_type == "GDD")
	{
	    GD_contribution << parameters("GD CONTRIBUTION");
	    merit_power << parameters("MERIT POWER");

	    target_EF << parameters("TARGET EF");
	    EF_tolerance << parameters("EF TOLERANCE");
	    EF_reserve << parameters("EF RESERVE");
	    dx << parameters("STEP");

	    target = new GDDTarget(number_of_frequencies,
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
				   number_of_bounces,
				   target_EF,
				   EF_tolerance,
				   EF_reserve,
				   dx,  // STEP
				   minimal_thickness,
				   0.0, /* individualism */
				   merit_power,
				   adaptive_dispersion);
	}
	else if (target_type == "pulse")
	{
	    str << parameters("PROBE PULSE");
	    if (! str.size())
	    {
		cerr << "\
ERROR: the target type is set to 'pulse', but the probe pulse is not\n\
defined. Assign a valid value to the option 'PROBE PULSE'.";
		return 1;
	    }
	    target = new PulseTarget(number_of_frequencies,
				     omega_min, omega_max,
				     variables.ProbePulseSpectrum(),
				     variables.ProbePulseSpectralPhase(),
				     variables.TargetDispersion(),
				     0.2, // time resolution
				     number_of_bounces,
				     minimal_thickness,
				     0.0, /* individualism */
				     adaptive_dispersion);
	}
	else
	{
	    cerr << "unknown optimisation target: " << target_type << "\n"
		 << "correct options: GDD, pulse." << endl;
	    return 1;
	}

	ConstructDesigns(designs, parameters, variables,
			 material_repository, *target);
	if (target) delete target;
	RenameDesignsToCWD(designs);
	// save the designs
	number_of_designs = designs.size();
	for (i=0; i<number_of_designs; i++) designs[i].SaveToFile();
    }
    catch(const string& error_message)
    {
	cerr << error_message << endl;
	return 1;
    }
    catch(const char* error_message)
    {
	cerr << error_message << endl;
	return 1;
    }
    catch(...)
    {
	cerr << "unknown exception" << endl;
	return 1;
    }

    return 0;
}

//-------------------------------------------------------------------------

void DesignNames(int argc, const char** argv, vector<string>& names)
{
    int i;
    string option;
    names.clear();
    for (i=0; i<argc; i++)
    {
	option = argv[i];
	if (option[0]=='-') continue;
	names.push_back(option);
    }
}

//-------------------------------------------------------------------------

void InitialiseParameters(Parameters& p)
{
    p.AddParameter("STEP", "10", "the step with which the stack has \
to be scanned (in nanometres)");
    p.AddParameter("TARGET EF", "0.7", "the target value of the \
	    normalised electric field");
    p.AddParameter("EF TOLERANCE", "0.05", "a tolerable average \
	    difference between the target and obtained EF");
    p.AddParameter("EF RESERVE", "0", "if at some particular wavelength \
	    the EF differs from the target EF by less than this value, \
	    then this point doesn't contribute to the merit function \
	    at all");
    p.AddParameter("VERBOSE", "no", "with this option the program \
prints some extra messages");
    p.AddParameter("INCIDENCE MEDIUM", "AIR",
		   "the name of the incidence medium");
    p.AddParameter("EXIT MEDIUM", "FUSI", "the name of the exit medium");
    p.AddParameter("ANGLE OF INCIDENCE", "0", "the angle of incidence");
    p.AddParameter("POLARISATION", "TM", "the polarisation of the beam (this \
parameter may have the following values: TM, p, TE, s, where 'TM' is the \
same as 'p' and 'TE' is the same as 's')");
    p.AddParameter("RESOLUTION", "100", "the number of frequencies at which \
the coating properties are evaluated during the optimisation stage \
(see also ANALYSIS RESOLUTION)");
    p.AddParameter("MINIMAL THICKNESS", "0",
		   "optimisation will try to avoid layers thinner than");
    p.AddParameter("MINIMAL WAVELENGTH", "500",
		   "the minimal wavelength of the \
spectral region used for optimisation and analysis (in nanometres)");
    p.AddParameter("MAXIMAL WAVELENGTH", "1000", 
		   "the maximal wavelength of the \
spectral region used for optimisation and analysis (in nanometres)");
    p.AddParameter("CENTRAL WAVELENGTH", "800",
		   "the wavelength at which dispersion is expanded into \
the Taylor series (in nanometres)");
    p.AddParameter("PROBE PULSE", "",
		   "the kind of probe pulse to be used for optimisation \
and analysis; if no value is specified for this parameter, the simulation \
of a probe pulse reflection is disabled; the probe pulse can be either \
'gaussian' or 'supergaussian' or a file name can be given in quotes. \
The file must contain the spectrum of the probe pulse: the wavelength in the \
first column, the spectral intensity in the second, and the optional \
spectral phase in radians in the third column. If the probe pulse is \
super-Gaussian, the word 'supergaussian' must be followed by a whitespace \
and an integer number, which is the order of the pulse ('1' corresponds to \
the usual Gaussian pulse, the higher the order, the more rectangular is \
the pulse); if the probe pulse is Gaussian or super-Gaussian, \
'PULSE CENTRAL WAVELENGTH' and 'PROBE PULSE FWHM' can specify its parameters");
    p.AddParameter("PROBE PULSE FWHM", "",
		   "the spectral width (FWHM) of the probe pulse \
(in nanometres); if the parameter is not initialised, a reasonable \
value will be used for the FWHM");
    p.AddParameter("PULSE CENTRAL WAVELENGTH", "", 
		   "the central wavelength of the probe pulse \
(in nanometres); if the parameter is not initialised, the value of \
CENTRAL WAVELENGTH will be used");
    p.AddParameter("MERIT POWER", "2", "the power of the merit function");
    p.AddParameter("TARGET DISPERSION", "0.0 0.0 0.0",
		   "this parameter specifies the target dispersion; it \
must be either a list of whitespace-separated numbers or a file name \
in quotes. In the first case the list of whitespace-separated numbers \
is the Taylor expansion of the target GDD (in fs^2) at the central \
wavelength: the first number is the GDD, the second number is the TOD \
and so on; if a file name in quotes is given, the file must contain the \
target GDD [fs^2] against the wavelength [nm]");
    p.AddParameter("ADAPTIVE DISPERSION", "no", 
		   "if the \"adaptive dispersion\" is activated, then each \
time a dispersion-based merit function is evaluated, the target dispersion \
will be adapted to the obtained one by multiplication by a proper factor");
    p.AddParameter("GD TOLERANCE", "3.0",
		   "a tolerable average difference between the target \
and obtained GD [fs]; the parameter must be either a number or a file name \
in quotes");
    p.AddParameter("GD RESERVE", "0", "if at some particular wavelength \
the GD differs from the target GD by less than this value, then this \
point doesn't contribute to the merit function at all; the parameter \
must be either a number or a file name in quotes");
    p.AddParameter("GDD TOLERANCE", "20",
		   "a tolerable average difference between the target \
and obtained GDD [fs^2]; the parameter must be either a number or a file name \
in quotes");
    p.AddParameter("GDD RESERVE", "0",
		   "if at some particular wavelength the GDD differs from \
the target GDD by less than this value, then this point doesn't contribute \
to the merit function at all; the parameter must be either a number or a \
file name in quotes");
    p.AddParameter("GD CONTRIBUTION", "0",
		   "a value between 0 and 1, which determines the weight \
of the GD against the GDD in the merit function");
    p.AddParameter("TARGET REFLECTANCE", "1.0",
		   "a value between 0 and 1.0 meaning the target \
reflectance (1 corresponds to the 100\% reflectance and 0 means full \
transmittance)");
    p.AddParameter("REFLECTANCE TOLERANCE", "0.01",
		   "a tolerable difference between the target value of \
the reflectance and the reflectance of the optimised design; the parameter \
must be either a number or a file name in quotes");
    p.AddParameter("REFLECTANCE RESERVE", "0",
		   "if at some particular wavelength the reflectance \
differs from the target reflectance by less than this value, then this \
point doesn't contribute to the merit function at all; the parameter must \
be either a number or a file name in quotes");
    p.AddParameter("NUMBER OF BOUNCES", "1", 
		   "the number of bounces off the mirror");
    p.AddParameter("LOCAL OPTIMISATION ALGORITHM", "gradient",
		   "the algorithm that will be used for the local \
optimisation (the valid names are 'simplex' and 'gradient')");
    p.AddParameter("TARGET TYPE", "GDD",
		   "the way how the merit function will be calculated \
(the valid options are 'GDD' and 'pulse')");
    p.AddParameter("SIMPLEX STEP", "10", "if the simplex algorithm is used \
for the local optimisation, this parameter specifies the variation \
of the layer thickness [nm] used to build the initial simplex");
    p.AddParameter("PATH TO MATERIALS", "",
		   "the data for each material is first searched in the \
current directory; if it's not found there, it is searched under this path");
    // parameters specific to 'acdinit'
    p.AddParameter("PARAMETRISATION ORDER", "8",
		   "the maximal number of parameters, which specify the \
thicknesses of layers with either low-or high refractive index");
}

//-------------------------------------------------------------------------

void InitialiseUnusedParameters(Parameters& p)
{
    p.AddParameter("COMPLEMENTARY", "", "");
    p.AddParameter("SENSITIVITY", "", "");
    p.AddParameter("ANALYSIS RESOLUTION", "", "");
    p.AddParameter("POPULATION SIZE", "", "");
    p.AddParameter("CONVERGENCE", "", "");
    p.AddParameter("NUMBER OF GENERATIONS", "", "");
    p.AddParameter("NUMBER OF RESULTS", "", "");
    p.AddParameter("INDIVIDUALISM", "", "");
    p.AddParameter("TECHNOLOGICAL ERROR", "", "");
    p.AddParameter("NUMBER OF TRIALS", "", "");
    p.AddParameter("MUTATION PROBABILITY", "", "");
    p.AddParameter("SAVE TRANSMITTANCE", "", "");
    p.AddParameter("SAVE ABSORPTION", "", "");
    p.AddParameter("SAVE PHASE", "", "");
}

//-------------------------------------------------------------------------

AbstractParameterConverter::AbstractParameterConverter()
{
    n_layers = 0;
    thicknesses = NULL;
}

AbstractParameterConverter::~AbstractParameterConverter()
{
    if (thicknesses) delete[] thicknesses;
}

void AbstractParameterConverter::
UpdateCoating(int number_of_parameters, Real* parameters, Coating& coating)
{
    int new_number_of_layers = coating.NumberOfVariableLayers();
    if (! thicknesses) thicknesses = new Real[new_number_of_layers];
    else if (n_layers < new_number_of_layers)
    {
	delete[] thicknesses;
	thicknesses = new Real[new_number_of_layers];
    }
    n_layers = new_number_of_layers;
    LayerThicknesses(number_of_parameters, parameters, n_layers, thicknesses);
    coating.SetVariableThicknesses(thicknesses);
}


//-------------------------------------------------------------------------

ChebyshevParameterConverter::ChebyshevParameterConverter() :
    AbstractParameterConverter() {}

void ChebyshevParameterConverter::
LayerThicknesses(int number_of_parameters, Real* parameters,
		 int number_of_layers, Real *layer_thicknesses)
{   /* The vector of parameters is interpreted as weight factors of
       Chebyshev polynomials of increasing orders, while even elements
       of 'parameters' parametrise layers with even indexes, and odd
       elements parametrise layers with odd indexes. */
    int i, j, n;
    Real x, thickness;

    for (i=0; i<number_of_layers; i++)
    {
	// the Chebyshev polynomials are defined in [-1:1]
	x = 2*i/Real(number_of_layers-1) - 1;
	thickness = 0;
	for (j=i%2, n=0; j<number_of_parameters; j+=2, n++)
	    thickness += parameters[j]*cos(n * acos(x));
	layer_thicknesses[i] = abs(thickness);
    }
    // remark: I could do it faster using the Clenshaw recursion
}

//-------------------------------------------------------------------------

void ConstructDesigns(vector<Design>& designs, Parameters& parameters,
		      Variables& variables,
		      MaterialRepository& material_repository,
		      AbstractTarget& target)
{
    const Real gradient_step = 1e-5;
    const Real simplex_step = 0.1;
    Real* optimisation_parameters;
    int i, n, n_parameters_max;
    int N = designs.size();
    bool verbose;
    Real lambda_min, lambda_max, helper;
    Real n_even_lambda_min, n_odd_lambda_min;
    Real n_even_lambda_max, n_odd_lambda_max;
    Real d1_even, d1_odd, d2_even, d2_odd;
    string material_name;
    enum Algorithm {gradient, simplex};
    Algorithm algorithm;
    string str;
    ChebyshevParameterConverter parameter_converter;
    AbstractDispersion& target_dispersion = variables.TargetDispersion();

    n_parameters_max << parameters("PARAMETRISATION ORDER");
    if (n_parameters_max < 2)
    {
	cerr << "\
WARNING: The parameter 'PARAMETRISATION ORDER' must be >= 2." << endl;
	n_parameters_max = 2;
    }
    optimisation_parameters = new Real[2*(n_parameters_max+1)];

    verbose << parameters("VERBOSE");
    str << parameters("LOCAL OPTIMISATION ALGORITHM");
    if (str=="gradient") algorithm = gradient;
    else if (str=="simplex") algorithm = simplex;
    else
    {
	cerr << "unknown local optimisation algorithm: " << str << endl;
	exit(1);
    }

    Coating coating;
    Coating::Parameters coating_parameters;
    coating_parameters.polarisation << parameters("POLARISATION");
    coating_parameters.angle_of_incidence << parameters("ANGLE OF INCIDENCE");
    coating_parameters.angle_of_incidence *= M_PI/180.0;
    coating_parameters.incidence_medium_name << parameters("INCIDENCE MEDIUM");
    coating_parameters.exit_medium_name << parameters("EXIT MEDIUM");
    coating_parameters.N << parameters("RESOLUTION");
    lambda_min << parameters("MINIMAL WAVELENGTH");
    lambda_max << parameters("MAXIMAL WAVELENGTH");
    coating_parameters.omega_min = 2.0*M_PI*SPEED_OF_LIGHT / lambda_max;
    coating_parameters.omega_max = 2.0*M_PI*SPEED_OF_LIGHT / lambda_min;
    coating.SetParameters(coating_parameters, material_repository);

    // if a positive-GDD mirror has to be designed, swap lambda_min and
    // lambda_max
    if (target_dispersion.GD(coating_parameters.omega_min) <
	target_dispersion.GD(coating_parameters.omega_max))
    {
	helper = lambda_min;
	lambda_min = lambda_max;
	lambda_max = helper;
    }

    // proceed the template designs one by one
    for (i=0; i<N; i++)
    {
	coating.ImportDesign(designs[i], material_repository);
	if (verbose) cout << "constructing " << designs[i].name << endl;
	// initialise the optimisation parameters
	material_name = designs[i].MaterialName(designs[i].NumberOfLayers()-1);
	n_even_lambda_min = real(
	    material_repository.RefractiveIndex(material_name, lambda_min));
	n_even_lambda_max = real(
	    material_repository.RefractiveIndex(material_name, lambda_max));
	material_name = designs[i].MaterialName(designs[i].NumberOfLayers()-2);
	n_odd_lambda_min = real(
	    material_repository.RefractiveIndex(material_name, lambda_min));
	n_odd_lambda_max = real(
	    material_repository.RefractiveIndex(material_name, lambda_max));
	d1_even = 0.25*lambda_min/n_even_lambda_min;
	d1_odd = 0.25*lambda_min/n_odd_lambda_min;
	d2_even = 0.25*lambda_max/n_even_lambda_max;
	d2_odd = 0.25*lambda_max/n_odd_lambda_max;
	optimisation_parameters[0] = 0.5*(d1_even + d2_even);
	optimisation_parameters[1] = 0.5*(d1_odd + d2_odd);
	optimisation_parameters[2] = 0.5*(d2_even - d1_even);
	optimisation_parameters[3] = 0.5*(d2_odd - d1_odd);
	for (n=4; n<=2*n_parameters_max; n++) optimisation_parameters[n] = 0;
	// optimise the design increasing the number of parameters
	for (n=2; n<=n_parameters_max; n++)
	{
	    if (verbose)
	    {
		cout << "\t" << n << " parameters: ";
		parameter_converter.UpdateCoating(2*n,
						  optimisation_parameters,
						  coating);
		cout << target.MeritFunction(coating) << " -> ";
		cout.flush();
	    }
	    switch(algorithm)
	    {
		case gradient:
		    ParametricCGoptimisation(2*n, optimisation_parameters,
					     parameter_converter,
					     coating, target, gradient_step);
		    break;
		case simplex:
		    ParametricSimplexOptimisation(2*n, optimisation_parameters,
						  parameter_converter,
						  coating, target,
						  simplex_step);
		    break;
	    }
	    if (verbose) cout << target.MeritFunction(coating) << endl;
	}
	coating.ExportDesign(designs[i]);
    }
    delete[] optimisation_parameters;
}

//-------------------------------------------------------------------------

double funk1(const gsl_vector *arguments, void *params)
{
    int i;
    if (optimised_coating && global_parameter_converter &&
	global_target && global_optimisation_parameters)
    {
	for (i=0; i< global_n_parameters; i++)
	    global_optimisation_parameters[i] = gsl_vector_get(arguments,i);
	global_parameter_converter->
	    UpdateCoating(global_n_parameters,
			  global_optimisation_parameters,
			  *optimised_coating);
	return global_target->MeritFunction(*optimised_coating);
    }
    else throw "ERROR in funk1";
    return 0.0;
}


void dfunk(const gsl_vector *arguments, void *params, gsl_vector *gradient)
{
    int i, n;
    Real original_value, v1, v2;
    if (optimised_coating && global_parameter_converter &&
	global_target && global_optimisation_parameters)
    {
	n = optimised_coating->NumberOfVariableLayers();
	for (i=0; i< global_n_parameters; i++)
	    global_optimisation_parameters[i] = gsl_vector_get(arguments, i);
	for (i=0; i< global_n_parameters; i++)
	{
	    // calculate the gradient with respect to the i-th parameter
	    original_value = global_optimisation_parameters[i];
	    global_optimisation_parameters[i] -= global_gradient_step;
	    global_parameter_converter->
		UpdateCoating(global_n_parameters,
			      global_optimisation_parameters,
			      *optimised_coating);
	    v1 = global_target->MeritFunction(*optimised_coating);
	    global_optimisation_parameters[i] = original_value +
		global_gradient_step;
	    global_parameter_converter->
		UpdateCoating(global_n_parameters,
			      global_optimisation_parameters,
			      *optimised_coating);
	    v2 = global_target->MeritFunction(*optimised_coating);
	    global_optimisation_parameters[i] = original_value;
	    gsl_vector_set(gradient, i, (v2 - v1) / (2.0*global_gradient_step));
	}
    }
    else throw "ERROR in 'dfunk'";
}

void fdfunk(const gsl_vector *arguments, void *params, double *f, gsl_vector *gradient)
{
    *f = funk1(arguments, params);
    dfunk(arguments, params, gradient);
}

//-------------------------------------------------------------------------


void ParametricCGoptimisation(int n, Real* parameters,
			      AbstractParameterConverter& converter,
			      Coating& coating, AbstractTarget& target,
			      Real step)
{
    const int params = NULL;    // GSL requires such parameters. Not needed
   				// in our program, though.

    const gsl_multimin_fdfminimizer_type *T; 
    gsl_multimin_fdfminimizer *gradient_minimizer = NULL;
    gsl_vector *x;
    gsl_multimin_function_fdf merit_min;

    int status;
    int i;
    const float ftol = 1e-4; // the fractional tolerance

    float* p; // the coordinates of a point in n-dimensional space
    size_t iter = 0; // the number of iterations that were performed

    float fret; // the minimum value of the function
    if (n<=0) return;

    // initialise the global variables
    global_target = &target;
    optimised_coating = &coating;
    global_parameter_converter = &converter;
    global_n_parameters = n;
    global_optimisation_parameters = new Real[n];
    global_gradient_step = step;

    // initialise the starting point
    p = new float[n+1];

    merit_min.n = n;
    merit_min.f = funk1;
    merit_min.df = dfunk;
    merit_min.fdf = fdfunk;
    merit_min.params = params;

    x = gsl_vector_alloc(n);
    for (i=0; i<n; i++) gsl_vector_set(x, i, parameters[i]);

    T = gsl_multimin_fdfminimizer_conjugate_fr;
    gradient_minimizer = gsl_multimin_fdfminimizer_alloc(T, n);

    // optimise
    gsl_multimin_fdfminimizer_set(gradient_minimizer, &merit_min, x, step, ftol);
    
    do {
	iter++;
	status = gsl_multimin_fdfminimizer_iterate(gradient_minimizer);
	
	if (status)
	    break;
	
	status = gsl_multimin_test_gradient(gradient_minimizer->gradient, ftol);
	
	/*if (status == GSL_SUCCESS)
	{
	    printf ("Minimum found.\n");
	}*/
    } while (status == GSL_CONTINUE && iter < 10000);
   
    for (i=0; i<n; i++)
	parameters[i] = gsl_vector_get(gradient_minimizer->x, i);
    
    fret = funk1(gradient_minimizer->x, params); // to ensure that 'coating' is updated
    
    // clean up
    delete[] p;
    delete[] global_optimisation_parameters;
    global_target = NULL;
    optimised_coating = NULL;
    global_parameter_converter = NULL;
    global_optimisation_parameters = NULL;
    global_gradient_step = 0;
}

//-------------------------------------------------------------------------

void ParametricSimplexOptimisation(int n, Real* parameters,
				   AbstractParameterConverter& converter,
				   Coating& coating, AbstractTarget& target,
				   Real step)
{
    const int params = NULL;    // GSL requires such parameters. Not needed
   				// in our program, though.

    const gsl_multimin_fminimizer_type *T = 
	gsl_multimin_fminimizer_nmsimplex2;
    gsl_multimin_fminimizer *simplex = NULL;
    gsl_vector *steps, *x;
    gsl_multimin_function merit_min;

    size_t iter = 0;
    int status;
    double size;

    const float ftol = 1e-4; // the fractional tolerance
    float fret;
    float* p; 
    int i;
    if (n<=0) return;

    // initialise the global variables
    global_target = &target;
    optimised_coating = &coating;
    global_parameter_converter = &converter;
    global_n_parameters = n;
    global_optimisation_parameters = new Real[n];

    p = new float[n];

    // initialise the simplex
    x = gsl_vector_alloc(n);
    for (i=0; i<n; i++) gsl_vector_set(x, i, parameters[i]);

    steps = gsl_vector_alloc(n);
    gsl_vector_set_all(steps, step);

    merit_min.n = n;
    merit_min.f = funk1;
    merit_min.params = params;

    simplex = gsl_multimin_fminimizer_alloc(T, n);
    gsl_multimin_fminimizer_set(simplex, &merit_min, x, steps);

    do {
	iter++;
	status = gsl_multimin_fminimizer_iterate(simplex);

	if (status) break;

	size = gsl_multimin_fminimizer_size(simplex);
	status = gsl_multimin_test_size(size, ftol);
	
	/*if (status == GSL_SUCCESS)
	{
	    printf("Converged to minumum.\n");
	}*/

    } while (status == GSL_CONTINUE && iter < 1000);
    
    for (i=0; i<n; i++)
	parameters[i] = gsl_vector_get(simplex->x, i);
    
    fret = funk1(simplex->x, params); // to ensure that 'coating' is updated

    // clean up
    delete[] p;
    delete[] global_optimisation_parameters;
    global_target = NULL;
    optimised_coating = NULL;
    global_parameter_converter = NULL;
    global_optimisation_parameters = NULL;
}

//-------------------------------------------------------------------------
