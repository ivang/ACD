// Calculates the electromagnetic field in the given layer stack

#include "common.hh"
#include "parameters.hh"
#include "design.hh"
#include "coating.hh"
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

//------------------------------------------------------------------------- 

void Usage(Parameters& parameters)
{
    cout << "USAGE: acdfield <options> <list of designs>\n";
    parameters.PrintHelp();
}

void DesignNames(int argc, const char** argv, vector<string>& names);
// extracts the design names from the command line

void InitialiseParameters(Parameters& parameters);

void InitialiseUnusedParameters(Parameters& parameters);
/* just to suppress the warnings */

void Analyse(vector<Design>& designs, Parameters& parameters,
	     MaterialRepository& material_repository);

int SaveField(const int Nx, const int Ny, const Real dx, 
	const Real lambda_min, const Real lambda_max, 
	Real** Z, const char* file_name);
/* returns 0 if the data have been successfully saved */

int Visualise(const int Nx, const int Ny, Real** Z,
	       const char* file_name, Real threshold=1e-8);
/* returns 0 if the visualisation has been successfully saved */

//-------------------------------------------------------------------------

int main(int argc, const char** argv)
{
    int i, number_of_designs;
    vector<string> design_names;
    string str, str1;

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
    try
    {
	str << parameters("ANALYSIS RESOLUTION");
	if (str.size() != 0)
	    parameters.SetParameterValue("RESOLUTION", str.c_str());
    }
    catch(...) {}

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
	RenameDesignsToCWD(designs);

	// analyse the designs
	Analyse(designs, parameters, material_repository);
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

void InitialiseParameters(Parameters& p)
{
    p.AddParameter("STEP", "10", "the step with which the stack has \
to be scanned (in nanometres)");
    p.AddParameter("THRESHOLD", "1e-4", "the lowest intensity, which \
appears in plots (the highest intensity being equal to 1)");
    p.AddParameter("SAVE DATA", "no", "the flag specifies whether \
the field analysis has to be saved in data files");
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
    p.AddParameter("MINIMAL WAVELENGTH", "500",
		   "the minimal wavelength of the \
spectral region used for optimisation and analysis (in nanometres)");
    p.AddParameter("MAXIMAL WAVELENGTH", "1000", 
		   "the maximal wavelength of the \
spectral region used for optimisation and analysis (in nanometres)");
    p.AddParameter("PATH TO MATERIALS", "",
		   "the data for each material is first searched in the \
current directory; if it's not found there, it is searched under this path");
}

//-------------------------------------------------------------------------

void InitialiseUnusedParameters(Parameters& p)
{
    p.AddParameter("ANALYSIS RESOLUTION", "", "");
    p.AddParameter("COMPLEMENTARY", "", "");
    p.AddParameter("NUMBER OF BOUNCES", "", "");
    p.AddParameter("CENTRAL WAVELENGTH", "", "");
    p.AddParameter("TARGET REFLECTANCE", "", "");
    p.AddParameter("REFLECTANCE TOLERANCE", "", "");
    p.AddParameter("REFLECTANCE RESERVE", "", "");
    p.AddParameter("TARGET DISPERSION", "", "");
    p.AddParameter("ADAPTIVE DISPERSION", "", "");
    p.AddParameter("GD CONTRIBUTION", "", "");
    p.AddParameter("GD TOLERANCE", "", "");
    p.AddParameter("GDD TOLERANCE", "", "");
    p.AddParameter("INDIVIDUALISM", "", "");
    p.AddParameter("MINIMAL THICKNESS", "", "");
    p.AddParameter("TARGET TYPE", "", "");
    p.AddParameter("POPULATION SIZE", "", "");
    p.AddParameter("NUMBER OF GENERATIONS", "", "");
    p.AddParameter("MUTATION PROBABILITY", "", "");
    p.AddParameter("NUMBER OF RESULTS", "", "");
    p.AddParameter("SENSITIVITY", "", "");
    p.AddParameter("NUMBER OF TRIALS", "", "");
    p.AddParameter("TECHNOLOGICAL ERROR", "", "");
    p.AddParameter("PROBE PULSE", "", "");
    p.AddParameter("SAVE TRANSMITTANCE", "", "");
    p.AddParameter("SAVE ABSORPTION", "", "");
    p.AddParameter("SAVE PHASE", "", "");
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

void Analyse(vector<Design>& designs, Parameters& parameters,
	     MaterialRepository& material_repository)
{
    int i, k, m, N, Nx;
    int n = designs.size();
    bool verbose, save_data;
    Real dx, lambda_min, lambda_max, threshold;
    Complex *E, *H, *eta;
    Real **Z;
    string str;

    // initialise the coating
    Coating coating;
    Coating::Parameters coating_parameters;
    try
    {
	verbose << parameters("VERBOSE");
	save_data << parameters("SAVE DATA");
	dx << parameters("STEP");
	threshold << parameters("THRESHOLD");
	coating_parameters.polarisation << parameters("POLARISATION");
	coating_parameters.angle_of_incidence <<
	    parameters("ANGLE OF INCIDENCE");
	coating_parameters.angle_of_incidence *= M_PI/180.0;
	coating_parameters.incidence_medium_name <<
	    parameters("INCIDENCE MEDIUM");
	coating_parameters.exit_medium_name << parameters("EXIT MEDIUM");
	coating_parameters.N << parameters("RESOLUTION");
	lambda_min << parameters("MINIMAL WAVELENGTH");
	lambda_max << parameters("MAXIMAL WAVELENGTH");
    }
    catch (const string& s)
    {
	cerr << "FAILED to get one of the important parameters:\n"
	     << s << endl;
	return;
    }
    if (dx<Real(0))
    {
	cerr << "ERROR: 'STEP' must be positive" << endl;
	return;
    }
    coating_parameters.omega_min = 2.0*M_PI*SPEED_OF_LIGHT / lambda_max;
    coating_parameters.omega_max = 2.0*M_PI*SPEED_OF_LIGHT / lambda_min;
    coating.SetParameters(coating_parameters, material_repository);

    // calculate the electric-field intensity
    N = coating_parameters.N;
    E = new Complex[N];
    H = new Complex[N];
    eta = new Complex[N];
    for (m=0; m<n; m++)
    {
	coating.ImportDesign(designs[m], material_repository);
	Nx = int(ceil(coating.StackThickness()/dx));
	Z = new Real*[Nx];
	coating.EFieldIntensity(N, Nx, dx, Z);

	if (save_data)
	{
	    str = designs[m].name + "_efi.dat";
	    k = SaveField(Nx, N, dx, lambda_min, lambda_max, Z, str.c_str());
	    if (verbose && k==0)
		cout << str << " saved" << endl;
	}
	str = designs[m].name + "_efi.ppm";
	k = Visualise(Nx, N, Z, str.c_str(), threshold);
	if (verbose && k==0)
	    cout << str << " saved" << endl;
	for (i=0; i<Nx; i++) delete[] Z[i];
	delete[] Z;
    } // for (m=0; m<n; m++)

    delete[] E;
    delete[] H;
    delete[] eta;
}

//-------------------------------------------------------------------------

int SaveField(const int Nx, const int Ny, const Real dx, 
	const Real lambda_min, const Real lambda_max, 
	Real** Z, const char* file_name)
{
    int i, j;
    FILE* file = fopen(file_name, "w");
    const Real dlambda = (lambda_max - lambda_min)/Ny;

    if (! file)
    {
	cerr << "FAILED to open " << file_name << endl;
	return -1;
    }
    for (i=0; i<Nx; i++)
    {
	for (j=0; j<Ny; j++) 
	{
	    fprintf(file, "%5.1f %5.1f %.4e\n",
		    i * dx, lambda_min + j * dlambda, Z[i][j]);
	}
	fprintf(file, "\n");
    }
    fclose(file);
    return 0;
}

//-------------------------------------------------------------------------

// r,g,b values are from 0 to 1
// h = [0,360], s = [0,1], v = [0,1]
//              if s == 0, then h = -1 (undefined)
void RGBtoHSV(float r, float g, float b, float *h, float *s, float *v);
void HSVtoRGB(float *r, float *g, float *b, float h, float s, float v);

int Visualise(const int Nx, const int Ny, Real** Z,
	       const char* file_name, Real threshold)
{
    int i, j;
    Real z, z_max;
    char pixel[3];
    float r, g, b, h, s, v;
    Real log_threshold = log(threshold);
    FILE* file;

    // find the maximum of Z
    z_max = 0;
    for (i=0; i<Nx; i++)
    {
	for (j=0; j<Ny; j++)
	{
	    if (z_max < Z[i][j]) z_max = Z[i][j];
	}
    }

    // make a plot
    file = fopen(file_name, "w");
    if (! file)
    {
	cerr << "FAILED to open " << file_name << endl;
	return -1;
    }

    fprintf(file, "%s\n%d %d\n%d\n", "P6", Nx, Ny, 255);
    for (j=Ny-1; j>=0; j--)
    {
	for (i=0; i<Nx; i++)
	{
	    z = Z[i][j] / z_max;
	    if (z <= threshold)
	    {
		r = 0.0;
		g = 0.0;
		b = 0.0;
	    }
	    else
	    {
		v = 1-log(z)/log_threshold;
		h = 300.0 * j / Real(Ny-1);
		s = 1;
		HSVtoRGB(&r, &g, &b, h, s, v);
	    }
	    pixel[0] = int(round(r*255));
	    pixel[1] = int(round(g*255));
	    pixel[2] = int(round(b*255));
	    fwrite(pixel, sizeof(pixel[1]), 3, file);
	}
    }
    fclose(file);
    return 0;
}

//-------------------------------------------------------------------------

void RGBtoHSV( float r, float g, float b, float *h, float *s, float *v )
{
  float min, max, delta;
  if (r > g) min = g;
  else min = r;
  if (min > b) min = b;
  if (r > g) max = r;
  else max = g;
  if (max < b) max = b;
  *v = max;                               // v
  delta = max - min;
  if( max != 0 )
    *s = delta / max;               // s
  else {
    // r = g = b = 0                // s = 0, v is undefined
    *s = 0;
    *h = -1;
    return;
  }
  if( r == max )
    *h = ( g - b ) / delta;         // between yellow & magenta
  else if( g == max )
    *h = 2 + ( b - r ) / delta;     // between cyan & yellow
  else
    *h = 4 + ( r - g ) / delta;     // between magenta & cyan
  *h *= 60;                               // degrees
  if( *h < 0 )
    *h += 360;
}

void HSVtoRGB( float *r, float *g, float *b, float h, float s, float v )
{
  int i;
  float f, p, q, t;
  if( s == 0 ) {
    // achromatic (grey)
    *r = *g = *b = v;
    return;
  }
  h /= 60;                        // sector 0 to 5
  i = int(floor(h));
  f = h - i;                      // factorial part of h
  p = v * ( 1 - s );
  q = v * ( 1 - s * f );
  t = v * ( 1 - s * ( 1 - f ) );
  switch( i ) {
  case 0:
    *r = v;
    *g = t;
    *b = p;
    break;
  case 1:
    *r = q;
    *g = v;
    *b = p;
    break;
  case 2:
    *r = p;
    *g = v;
    *b = t;
    break;
  case 3:
    *r = p;
    *g = q;
    *b = v;
    break;
  case 4:
    *r = t;
    *g = p;
    *b = v;
    break;
  default:                // case 5:
    *r = v;
    *g = p;
    *b = q;
    break;
  }
}
