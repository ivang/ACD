#include "common.hh"
#include "parameters.hh"
#include "material.hh"
#include "design.hh"
#include "coating.hh"
#include "optimisation.hh"
#include "analysis.hh"
#include "random.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <string.h>

//-------------------------------------------------------------------------

const Real mutation_std = 20.0; // nm

//-------------------------------------------------------------------------

void Usage(Parameters& parameters)
{
    string s;
    cout << "Advanced Coating Designer\n\n\
USAGE: acd [command] [options] <designs>\n" << endl;
    s = "where the 'command' has to be either 'l' (local optimisation) or \
'g' (global optimisation). If no command is given, the designs are only \
analysed without optimisation. If the global optimisation is requested, \
the designs are added to the initial population. If the local optimisation \
is requested and the option COMPLEMENTARY is set to 'no' or 'false', then \
all the designs are locally optimised one by one. If the option is set \
to 'yes' or 'true', all the designs are assumed to build a system, the \
properties of which can be optimised and/or analysed. The global \
optimisation of a set of complementary mirrors is not implemented.";
    PrintFormatted(cout, s, 75, "");
    cout << "  ";
    parameters.PrintHelp();
}

void ComplementaryOptimisation(vector<Design>& designs, Parameters& parameters,
			  MaterialRepository& material_repository,
			  AbstractTarget& target);

void LocalOptimisation(vector<Design>& designs, Parameters& parameters,
		       MaterialRepository& material_repository,
		       AbstractTarget& target);

void GlobalOptimisation(vector<Design>& designs, Parameters& parameters, 
			MaterialRepository& material_repository,
			AbstractTarget& target);

void DesignNames(int argc, const char** argv, vector<string>& names);
// extracts the design names from the command line

void InitialiseParameters(Parameters& parameters);

//-------------------------------------------------------------------------

int main(const int argc, const char** argv)
{
    int i, number_of_designs, analysis_resolution;
    vector<string> design_names;
    string str, str1;
    AbstractTarget* target = NULL;
    int number_of_frequencies, number_of_bounces, merit_power;
    Real lambda0, lambda_min, lambda_max, omega0, omega_min, omega_max;
    Real minimal_thickness, individualism;
    bool adaptive_dispersion, complementary_coatings;
    Real GD_contribution, time_resolution;
    Real target_EF, EF_tolerance, EF_reserve, dx;

    enum OptimisationMode {no_optimisation, local, global};
    OptimisationMode optimisation_mode;

    InitialiseRandomNumberGenerator();

    // initalise the parameters
    Parameters parameters;
    InitialiseParameters(parameters);
    if (argc<2) Usage(parameters);
    str = argv[1];
    if (str=="-h" || str=="--help") Usage(parameters);

    // how do we have to optimise the designs?
    optimisation_mode = no_optimisation;
    if (strlen(argv[1]) == 1)
    {
	if (argv[1][0]=='l') optimisation_mode = local;
	else if (argv[1][0]=='g') optimisation_mode = global;
    }

    // read the parameter file
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
	if (optimisation_mode == no_optimisation)
	    parameters.ParseCommandLine(argc-1, &argv[1]);
	else parameters.ParseCommandLine(argc-2, &argv[2]);
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
	complementary_coatings << parameters("COMPLEMENTARY");
	if (complementary_coatings)
	    individualism << parameters("INDIVIDUALISM");
    }
    catch (const string& s)
    {
	cerr << "FAILED to get one of the important parameters:\n"
	     << s << endl;
	return 1;
    }
    if (complementary_coatings && optimisation_mode==global)
    {
	cerr << "\
WARNING: the option COMPLEMENTARY does not have any effect on the global\n\
optimisation";
	complementary_coatings = false;
    }

    // there are a few parameters, which are allowed to be uninitialised
    // at this point; we must initialise them now
    str << parameters("ANALYSIS RESOLUTION");
    if (str.size() == 0)
    {
	str << parameters("RESOLUTION");
	parameters.SetParameterValue("ANALYSIS RESOLUTION", str.c_str());
    }
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
    if (optimisation_mode == no_optimisation)
	DesignNames(argc-1, &argv[1], design_names);
    else DesignNames(argc-2, &argv[2], design_names);

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
	if (number_of_designs==1) complementary_coatings = false;

	if (optimisation_mode != no_optimisation)
	{
	    // create a target
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
				       dx,
				       minimal_thickness,
				       individualism,
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
		time_resolution = 0.2 / (omega_max - omega_min);
		target = new PulseTarget(number_of_frequencies,
					 omega_min, omega_max,
					 variables.ProbePulseSpectrum(),
					 variables.ProbePulseSpectralPhase(),
					 variables.TargetDispersion(),
					 time_resolution,
					 number_of_bounces,
					 minimal_thickness,
					 individualism,
					 adaptive_dispersion);
	    }
	    else
	    {
		cerr << "unknown optimisation target: " << target_type << "\n"
		     << "correct options: GDD, pulse." << endl;
		return 1;
	    }
	}

	// optimisation
	switch(optimisation_mode)
	{
	    case no_optimisation:
		break;
	    case local:
		if (complementary_coatings)
		{
		    ComplementaryOptimisation(designs, parameters,
					      material_repository, *target);
		}
		else
		{
		    LocalOptimisation(designs, parameters,
				      material_repository, *target);
		}
		break;
	    case global:
		GlobalOptimisation(designs, parameters, material_repository,
				   *target);
		break;
	    default:
		cerr << "ERROR: unknown optimisation mode" << endl;
		return 1;
	}
	if (target) delete target;
	RenameDesignsToCWD(designs);
	if (optimisation_mode != no_optimisation)
	{   // save the designs
	    number_of_designs = designs.size();
	    for (i=0; i<number_of_designs; i++) designs[i].SaveToFile();
	}
	// analyse the designs
	analysis_resolution << parameters("ANALYSIS RESOLUTION");
	Analyse(designs, complementary_coatings,
		analysis_resolution, material_repository,
		parameters, variables);
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

//=========================================================================

void LocalOptimisation(vector<Design>& designs, Parameters& parameters,
		       MaterialRepository& material_repository,
		       AbstractTarget& target)
{
    int i;
    int N = designs.size();
    bool verbose;
    Real lambda_min, lambda_max, simplex_step;
    enum Algorithm {gradient, simplex};
    Algorithm algorithm;
    string str;

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

    // optimise the designs one by one
    for (i=0; i<N; i++)
    {
	coating.ImportDesign(designs[i], material_repository);
	if (verbose)
	{
	    cout << "optimising " << designs[i].name << ": "
		 << target.MeritFunction(coating) << " -> ";
	    cout.flush();
	}
	switch(algorithm)
	{
	    case gradient:
		ConjugateGradientOptimisation(coating, target);
		break;
	    case simplex:
		simplex_step << parameters("SIMPLEX STEP");
		SimplexOptimisation(coating, target, simplex_step);
		break;
	}
	if (verbose) cout << target.MeritFunction(coating) << endl;
	coating.ExportDesign(designs[i]);
    }
}


void ComplementaryOptimisation(vector<Design>& designs, Parameters& parameters,
			       MaterialRepository& material_repository,
			       AbstractTarget& target)
{
    int i;
    int N = designs.size();
    bool verbose;
    Real lambda_min, lambda_max, simplex_step;
    enum Algorithm {gradient, simplex};
    Algorithm algorithm;
    string str;

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


    // create a vetor of coatings
    vector<Coating> coatings;
    for (i=0; i<N; i++)
    {
	coating.ImportDesign(designs[i], material_repository);
	coatings.push_back(coating);
    }
    if (verbose)
    {
	cout << "optimising the set of mirrors: "
	     << target.MeritFunction(coatings) << " -> ";
	cout.flush();
    }
    // optimise
    switch(algorithm)
    {
	case gradient:
	    ConjugateGradientOptimisation(coatings, target);
	    break;
	case simplex:
	    simplex_step << parameters("SIMPLEX STEP");
	    SimplexOptimisation(coatings, target, simplex_step);
	    break;
    }
    if (verbose) cout << target.MeritFunction(coatings) << endl;
    // export the optimised designs
    for (i=0; i<N; i++)
    {
	coatings[i].ExportDesign(designs[i]);
    }
}


void GlobalOptimisation(vector<Design>& designs, Parameters& parameters,
			MaterialRepository& material_repository,
			AbstractTarget& target)
{
    MemeticOptimiser optimiser;
    int i, n, number_of_generations;
    int number_of_designs = designs.size();
    bool verbose;
    verbose << parameters("VERBOSE");
    Real lambda_min, lambda_max, merit, convergence, mutation_probability;
    int population_size;
    population_size << parameters("POPULATION SIZE");
    number_of_generations << parameters("NUMBER OF GENERATIONS");
    convergence << parameters("CONVERGENCE");
    mutation_probability << parameters("MUTATION PROBABILITY");
    Design design;
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

    // check if all the designs are similar
    for (i=1; i<number_of_designs; i++)
    {
	if (! designs[0].isLike(designs[i]))
	{
	    cerr << "\
ERROR: in order to optimise designs globally, they all must have the same \n\
layer structure and identical layers with fixed thicknesses (if any). The\n\
designs '" << designs[0].name << "' and '" << designs[i].name
		 << "' are different in this respect" << endl;
	    exit(1);
	}
    }

    // create a populaton
    Population population(target, mutation_probability, mutation_std);
    for (i=0; i<number_of_designs; i++)
    {
	coating.ImportDesign(designs[i], material_repository);
        population.Add(coating);
    }
    if (number_of_designs<population_size)
    {
	if (verbose) cout << "generating random desings" << endl;
	for (i=number_of_designs; i < population_size; i++)
	    population.AddRandomDesign();
    }

    // if the global optimisation makes sense at all, do it
    if (coating.NumberOfVariableLayers() &&
	(number_of_generations>0 || convergence>Real(0)))
    {
	// improve the members of the populaton a bit
	if (verbose)
	    cout << "preliminary optimisation of the initial population:"
		 << endl;
	for (i=0; i<population.Size(); i++)
	{
	    if (verbose)
	    {
		cout << population.MeritFunction(i) << " --> ";
		cout.flush();
	    }
	    merit = population.RefineQuickly(i);
	    if (verbose) cout << merit << endl;
	}
	// global optimisation
	if (verbose)
	{
	    Real best_merit = population.BestMerit();
	    Real new_best_merit;
	    if (verbose) cout << "GLOBAL OPTIMISATION\n\
# column 1: generation;\n\
# column 2: the figure of merit of the best member.\n";
	    cout << setw(7) << 0 <<"    "
		 << setprecision(7) << best_merit << endl;
	    for (i=0; ; i++)
	    {
		optimiser.Evolve(population, 1, population_size);
		new_best_merit = population.BestMerit();
		if (best_merit > new_best_merit)
		{
		    best_merit = new_best_merit;
		    cout << setw(7) << i+1 <<"    "
			 << setprecision(7) << best_merit << endl;
		}
		if (number_of_generations>0 && i>=number_of_generations)
		    break;
		if (convergence>Real(0) &&
		    optimiser.ConvergenceFactor(population) <= convergence)
		    break;
	    }
	    cout << "the convergence factor of the last generation: "
		 << optimiser.ConvergenceFactor(population) << endl;
	}
	else
	{
	    if (convergence==Real(0) && number_of_generations>0)
	    {
		optimiser.Evolve(population, number_of_generations,
				 population_size);
	    }
	    else if (number_of_generations>0)
	    {
		for (i=0; i<number_of_generations; i++)
		{
		    optimiser.Evolve(population, 1, population_size);
		    if (optimiser.ConvergenceFactor(population) <= convergence)
			break;
		}
	    }
	    else if (convergence>Real(0))
	    {
		do
		{
		    optimiser.Evolve(population, 100, population_size);
		} while(optimiser.ConvergenceFactor(population) > convergence);
	    }
	} // if (verbose)
    } // if (coating.NumberOfVariableLayers() && ...
    // refine and export a few best members
    n << parameters("NUMBER OF RESULTS");
    population.Sort();
    population_size = population.Size();
    if (n > population_size) n = population_size;
    if (verbose) cout << "refining " << n << " best members... " << flush;
    optimiser.RefineBestMembers(population, n); // and maybe kill some of them
    if (verbose)
    {
	cout << "done" << endl;
	if (population_size != population.Size())
	    cout << population_size-population.Size()
		 << " members have been killed" << endl;
	if (n > population.Size())
	    cout << "only " << population.Size()
		 << " designs will be saved" << endl;
    }
    population_size = population.Size();
    if (n > population_size) n = population_size;
    population.Sort();
    designs.clear();
    ostringstream os;
    for (i=0; i<n; i++)
    {
        population.Get(i, coating);
	coating.ExportDesign(design);
	os.str("");
	os << "design";
	if (n>9) os << setw(2) << setfill('0');
	else if (n>99) os << setw(3) << setfill('0');
	else if (n>999) os << setw(4) << setfill('0');
	os << i+1;
	design.name = os.str();
	designs.push_back(design);
    }
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
    p.AddParameter("COMPLEMENTARY", "no", "with this option the designs \
passed on the command line are assumed to build a set of complementary \
coatings");
    p.AddParameter("SENSITIVITY", "no",
		   "this option determines whether the sensitivity \
analysis has to be performed");
    p.AddParameter("INCIDENCE MEDIUM", "AIR",
		   "the name of the incidence medium");
    p.AddParameter("EXIT MEDIUM", "FUSI", "the name of the exit medium");
    p.AddParameter("ANGLE OF INCIDENCE", "0", "the angle of incidence");
    p.AddParameter("POLARISATION", "TM", "the polarisation of the beam (this \
parameter may have the following values: TM, p, TE, s, where 'TM' is the \
same as 'p' and 'TE' is the same as 's')");
    p.AddParameter("ANALYSIS RESOLUTION", "",
		   "the number of frequencies at which the coating \
properties are evaluated during the analysis stage; if this parameter \
is not specified, the resolution is given by the parameter RESOLUTION");
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
    p.AddParameter("POPULATION SIZE", "100", "the population size to be \
used for the global optimisation");
    p.AddParameter("CONVERGENCE", "0.0",
		   "the value which indicates the \
convergence of an evolutionary algorithm (0.01 is a reasonable value for \
this parameter, see also NUMBER OF GENERATIONS)");
    p.AddParameter("NUMBER OF GENERATIONS", "1000", 
		   "the maximal number of generations \
an evolutionary algorithm can go through (if this value is equal to \
zero, it means that only the 'CONVERGENCE' will determine \
when to stop the iterations)");
    p.AddParameter("NUMBER OF RESULTS", "10",
		   "the number of distinct designs \
to be refined and saved after the global optimisation");
    p.AddParameter("NUMBER OF BOUNCES", "1", 
		   "the number of bounces off the mirror");
    p.AddParameter("INDIVIDUALISM", "0.5",
		   "a number between 0 and 1, which defines how important \
it is that each mirror in a set of complementary mirrors has nice \
properties alone: if the parameter is equal to 0, only the properties of \
the whole system are considered; if it is equal to 1, the coatings \
are optimised independently of each other");
    p.AddParameter("TECHNOLOGICAL ERROR", "1.0", "the mean value of random \
technological perturbations used for the sensitivity analysis (in \
nanometres)");
    p.AddParameter("NUMBER OF TRIALS", "1000",
		   "the number of trials used to estimate the effect of \
perturbations of layer thicknesses due to technological errors");
    p.AddParameter("LOCAL OPTIMISATION ALGORITHM", "gradient",
		   "the algorithm that will be used for the local \
optimisation (the valid names are 'simplex' and 'gradient')");
    p.AddParameter("TARGET TYPE", "GDD",
		   "the way how the merit function will be calculated \
(the valid options are 'GDD' and 'pulse')");
    p.AddParameter("MUTATION PROBABILITY", "0.02", "mutation plays an \
important role in the global optimisation; this parameter must have a value \
between 0 and 1, which determines the probability of each 'gene' \
(layer thickness) to be mutated (randomly altered)");
    p.AddParameter("SIMPLEX STEP", "10", "if the simplex algorithm is used \
for the local optimisation, this parameter specifies the variation \
of the layer thickness [nm] used to build the initial simplex");
    p.AddParameter("PATH TO MATERIALS", "",
		   "the data for each material is first searched in the \
current directory; if it's not found there, it is searched under this path");
    p.AddParameter("SAVE TRANSMITTANCE", "", "if this option is set to \
\"yes\" or \"true\", the transmittance of the designs is calculated and \
saved; if it is set to \"no\" of \"false\", the transmittance is not saved; \
if this option is not set, the transmittance is only saved if one of the \
materials has a non-zero absorption or gain");
    p.AddParameter("SAVE ABSORPTION", "", "if this option is set to \
\"yes\" or \"true\", the absorption of the designs is calculated and \
saved; if it is set to \"no\" of \"false\", the absorption is not saved; \
if this option is not set, the absorption is only saved if one of the \
materials has a non-zero absorption or gain");
    p.AddParameter("SAVE PHASE", "no", "if this option is set to \
\"yes\" or \"true\", the phase shift upon reflection is calculated and \
saved; if this option is set to \"no\" of \"false\", the phase is not saved");
}

//-------------------------------------------------------------------------
