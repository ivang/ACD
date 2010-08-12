#include "parameters.hh"
#include "reader.hh"
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cmath>

                  //--------------------------------

void PrintFormatted(ostream& os, const string& str, int L, const string prefix)
{   // sends the given string to the given stream, trying not
    // to put more than 'L' symbols on a line and writing the 'prefix'
    // in the beginning of each line
    int n, s_length, prefix_length;
    istringstream is(str);
    string s;
    bool insert_whitespace = false;

    prefix_length = prefix.length();
    n = 0; // the number of symbols on the current line
    while (is >> s)
    {
	s_length = s.length();
	if (n==0 && prefix_length>0)
	{
	    os << prefix;
	    n = prefix_length;
	    insert_whitespace = false;
	}
	if (s_length+1+prefix_length > L)
	{   // a very long word
	    os << s;
	    n = 0;
	    continue;
	}
	if (n+s_length+1 <= L)
	{   // this word can be written on the current line
	    if (insert_whitespace)
	    {
		os << " ";
		n++;
	    }
	    os << s;
	    n += s_length;
	    insert_whitespace = true;
	}
	else
	{   // this word has to be written on the next line
	    os << "\n" << prefix << s;
	    n = prefix_length + s_length;
	    insert_whitespace = true;
	}
    }
    os << "\n";
}

//===========================================================================


AbstractParameter::AbstractParameter() {}

AbstractParameter::~AbstractParameter() {}

ConstantParameter::ConstantParameter(Real value) : value(value) {}

ConstantParameter::~ConstantParameter() {}

ConstantParameter::ConstantParameter(const ConstantParameter& p)
{
    value = p.value;
}

ConstantParameter& ConstantParameter::operator=(const ConstantParameter& p)
{
    value = p.value;
    return *this;
}

Real ConstantParameter::operator()(Real x) const
{
    return value;
}

//-------------------------------------------------------------------------

VariableParameter::VariableParameter()
{
    N = 0;
}

VariableParameter::~VariableParameter() {}

VariableParameter::VariableParameter(const VariableParameter& p)
{
    N = p.N;
    data = p.data;
}

VariableParameter& VariableParameter::operator=(const VariableParameter& p)
{
    N = p.N;
    data = p.data;
    return *this;
}

bool VariableParameter::isInitialised()
{
    return (N>0);
}

void VariableParameter::TakeData(const vector< pair<Real,Real> >& _data)
{
    data = _data;
    sort(data.begin(), data.end());
    N = data.size();
}

void VariableParameter::TakeData(vector<Real> X, vector<Real> Y)
{
    int i;
    if (X.size() != Y.size()) throw "ERROR in \
VariableParameter::TakeData(vector<Real>, vector<Real>)";
    data.clear();
    N = X.size();
    for (i=0; i<N; i++) data.push_back(make_pair(X[i], Y[i]));
    sort(data.begin(), data.end());
}

int VariableParameter::TakeData(const string& str)
{
    int i, n;
    string s;
    Real x, y;

    data.clear();
    N = 0;
    n = str.length() - 1;
    while (n >= 0 && (str[n] == '\n' || str[n] == '\0' ||
		      str[n] == ' ' || str[n] == '\t')) n--;
    if ((str[0] != '"' || str[n] != '"') &&
	(str[0] != '\'' || str[n] != '\''))
    {	// this is supposed to be just a number
	istringstream is(str);
	is >> y;
	if (is.fail())
	{
	    cerr << "ERROR: can't convert " << str << " to an integer number";
	    return 1;
	}
	x = 0.0;
	data.push_back(make_pair(x, y));
	N = 1;
	return 0;
    }
    // if we got to here, we are supposed to read the data from a file
    s = str.substr(1,n-1);
    ifstream data_file(s.c_str());
    if (! data_file)
    {
	cerr << "ERROR: can't to open the file \"" << s << "\"" << endl;
	return -1;
    }
    vector< vector<Real> > table;
    ReadTable(data_file, table);
    data_file.close();
    n = table.size();
    if (n < 1)
    {
	cerr << "ERROR: no data have been read from '" << s.c_str()
	     << "'" << endl;
	exit(1);
    }
    for (i=0; i<n; i++)
    {
	if (table[i].size() < 2)
	{
	    cerr << "ERROR: can't parse '" << s << "'" << endl;
	    exit(1);
	}
	data.push_back(make_pair(table[i][0], table[i][1]));
    }
    N = data.size();
    sort(data.begin(), data.end());
    return 0;
}

Real VariableParameter::operator()(Real x) const
{
    Real h;
    int k, k1, k2;
    if (N==0) return 0.0;
    if (N==1) return data[0].second;
    if (x < data[0].first) return data[0].second;
    if (x > data[N-1].first) return data[N-1].second;
    k1 = 0;
    k2 = N-1;
    // we will find the right place in the table by means of bisection
    while (k2-k1 > 1)
    {
	k = (k2+k1) / 2;
	if (data[k].first > x) k2 = k;
	else k1 = k;
    }
    // k1 and k2 now bracket the input value of x.
    h = data[k2].first - data[k1].first;
    if (h == 0.) throw "ERROR in VariableParameter::operator(): \
the X's must be distinct";
    return data[k1].second + (data[k2].second - data[k1].second) *
	(x - data[k1].first) / h;
}

//=========================================================================

FlexibleDispersion::FlexibleDispersion() {}

FlexibleDispersion::FlexibleDispersion(const vector<Real>& wavelengths, 
				       const vector<Real>& _GDD)
{
    Set(wavelengths, _GDD);
}

FlexibleDispersion::~FlexibleDispersion()
{
}

void FlexibleDispersion::Set(const vector<Real>& wavelengths,
			     const vector<Real>& _GDD)
{
    int i;
    vector<Real> omega, _GD, _phase;
    vector< pair<Real,Real> > data;
    int N = wavelengths.size();
    if (wavelengths.size() != _GDD.size())
	throw "ERROR in FlexibleDispersion::Set";
    omega.resize(N);
    _GD.resize(N);
    _phase.resize(N);
    for (i=0; i<N; i++) omega[i] = 2.0*M_PI*SPEED_OF_LIGHT / wavelengths[i];
    // sort the frequencies
    for (i=0; i<N; i++) data.push_back(make_pair(omega[i], _GDD[i]));
    sort(data.begin(), data.end());
    for (i=0; i<N; i++) omega[i] = data[i].first;
    // initialise the dispersion
    GDD_parameter.TakeData(data);
    _GD[0] = 0;
    for (i=1; i<N; i++)
	_GD[i] = _GD[i-1] + (data[i-1].second + data[i].second) *
	    (omega[i] - omega[i-1]) / 2;
    GD_parameter.TakeData(omega, _GD);
    _phase[0] = 0;
    for (i=1; i<N; i++)
    {
	_phase[i] = _phase[i-1] -
	    (_GD[i-1] + _GD[i]) * (omega[i] - omega[i-1]) / 2;
    }
    phase_parameter.TakeData(omega, _phase);
}

void FlexibleDispersion::Set(const char* file_name)
{
    int i, n;
    vector<Real> X, Y;
    vector< vector<Real> > table;
    ifstream file(file_name);
    if (! file)
    {
	cerr << "ERROR in FlexibleDispersion::Set: can't open the file "
	     << file_name << endl;
	exit(1);
    }
    ReadTable(file, table);
    file.close();
    n = table.size();
    if (n < 1)
    {
	cerr << "ERROR: no data have been read from '" << file_name
	     << "'" << endl;
	exit(1);
    }
    for (i=0; i<n; i++)
    {
	if (table[i].size() < 2)
	{
	    cerr << "ERROR in FlexibleDispersion::Set: \
can't read the target GDD from '" << file_name << "'" << endl;
	    exit(1);
	}
	X.push_back(table[i][0]);
	Y.push_back(table[i][1]);
    }
    Set(X, Y);
}

Real FlexibleDispersion::Phase(Real omega) const
{
    return phase_parameter(omega);
}

Real FlexibleDispersion::GD(Real omega) const
{
    return GD_parameter(omega);
}

Real FlexibleDispersion::GDD(Real omega) const
{
    return GDD_parameter(omega);
}

//=========================================================================

Parameter::Parameter(const Parameter& p)
{
    value = p.value;
    description = p.description;
}

const Parameter& Parameter::operator=(const Parameter& p)
{
    value = p.value;
    description = p.description;
    return *this;
}


int operator << (int& x, const Parameter& p)
{
    istringstream is(p.value);
    is >> x;
    if (is.fail())
    {
	string s = "can't convert " + p.value + " to an integer number";
	throw s;
    }
    return x;
}

float operator << (float& x, const Parameter& p)
{
    istringstream is(p.value);
    is >> x;
    if (is.fail())
    {
	string s = "can't convert " + p.value + " to an float";
	throw s;
    }
    return x;
}

double operator << (double& x, const Parameter& p)
{
    istringstream is(p.value);
    is >> x;
    if (is.fail())
    {
	string s = "can't convert " + p.value + " to double";
	throw s;
    }
    return x;
}

bool operator << (bool& x, const Parameter& p)
{
    string value(p.value);
    transform(value.begin(), value.end(), value.begin(), (int(*)(int))toupper);
    if (value == "YES" || value == "TRUE") x = true;
    else if (value == "NO" || value == "FALSE") x = false;
    else
    {
	string s = "can't interpret '" + p.value + 
	    "' as a boolean value";
	throw s;
    }
    return x;
}

Tpolarisation operator << (Tpolarisation& x, const Parameter& p)
{
    string value(p.value);
    transform(value.begin(), value.end(), value.begin(), (int(*)(int))toupper);
    if (value == "TM" || value == "P") x = TM;
    else if (value == "TE" || value == "S") x = TE;
    else
    {
	string s = "can't interpret '" + p.value + 
	    "' as a polarisation kind";
	throw s;
    }
    return x;
}

const string& operator << (string& x, const Parameter& p)
{
    x = p.value;
    return x;
}

ostream& operator<< (ostream& o, const Parameter& p)
{
    return o << p.value;
}

//=========================================================================

Parameters::Parameters(map<string,Parameter>& input_map)
{
    map<string,Parameter>::iterator it;
    string parameter_name;
    Parameter p;

    for (it=input_map.begin(); it!=input_map.end(); ++it)
    {
	// convert all the parameter names to uppercase
	parameter_name = it->first;
	transform(parameter_name.begin(), parameter_name.end(),
		  parameter_name.begin(), (int(*)(int)) toupper);
	// add the parameter
	parameter_map.insert(make_pair(parameter_name, it->second));
    }
}

//-------------------------------------------------------------------------


void Parameters::AddParameter(const char* name, const char* value,
			      const char* description)
{
    Parameter *p;
    string parameter_name(name);
    map<string,Parameter>::iterator it;
    transform(parameter_name.begin(), parameter_name.end(),
	      parameter_name.begin(), (int(*)(int)) toupper);
    // if the parameter is already in the map, just update its fields
    it = parameter_map.find(parameter_name);
    if (it != parameter_map.end())
    {
	p = &(it->second);
	p->value = value;
	p->description = description;
    }
    else // add the parameter to the map
    {
	Parameter p;
	p.value = value;
	p.description = description;
	parameter_map.insert(make_pair(parameter_name, p));
    }
}

//-------------------------------------------------------------------------


void Parameters::ParseFile(istream& is)
{
    const int max_line_length = 1024;
    int i, j, line_counter, str_length;
    char line[max_line_length];
    string str, parameter_name;
    string::size_type idx;
    map<string,Parameter>::iterator it;
    Parameter p;

    line_counter = 0;
    while (is.getline(line, max_line_length))
    {
	line_counter++;
	// if the first non-blank character is '#', it's a comment string
	// we also should check whether it's an empty string
	for (i=0; i<max_line_length && line[i]!='\0' && line[i]!='\n' &&
		 (line[i]==' ' || line[i]=='\t'); i++) {}
	if (i==max_line_length) continue;
	if (line[i]=='\0' || line[i]=='\n' || line[i]=='#') continue;
	// try to interpret the string
	str = line; // and we know that the first non-blank character has
		    // index 'i'
	idx = str.find("=", i);
	if (idx == string::npos)
	{
	    ostringstream os;
	    os << "Can't parse line number " << line_counter
	       << " of the parameter file:\n" << line << endl;
	    string s = os.str();
	    throw s;
	}
	// find the last non-blank character before the '=' sign
	for (j=idx-1; j>i && (str[j]==' ' || str[j]=='\t'); j--) {}
	// extract the parameter name
	parameter_name = str.substr(i, j-i+1);
	transform(parameter_name.begin(), parameter_name.end(),
		  parameter_name.begin(), (int(*)(int)) toupper);
	// extract the parameter value
	str_length = str.length();
	for (i=idx+1; i<str_length && (str[i]==' ' || str[i]=='\t'); i++) {}
	if (i==str_length) continue; // the parameter has no value
	for (j=i; j<str_length && str[j]!='#'; j++) {}
	for (j--; j>i && (str[j]==' ' || str[j]=='\t'); j--) {}
	p.value = str.substr(i, j-i+1);
	// do we know a parameter with such a name?
	it = parameter_map.find(parameter_name);
	if (it != parameter_map.end())
	{ // just update the value
	    (it->second).value = p.value;
	}
	else
	{ // issue a warning and add the parameter to the map
	    cerr << "WARNING: the parameter '" << parameter_name
		 << "' is not known." << endl;
	    p.description = "";
	    parameter_map.insert(make_pair(parameter_name, p));
	}
    } // while (is.getline(line, max_line_length))
}

//-------------------------------------------------------------------------


void Parameters::ParseCommandLine(int argc, const char** argv)
{
    int i, j, n;
    string option, parameter_name;
    string::size_type idx;
    map<string,Parameter>::iterator it;
    Parameter p;

    for (i=0; i<argc; i++)
    {
	option = argv[i];
	if (option=="-h" || option=="--help")
	{
	    PrintHelp();
	    exit(0);
	}
	if (option.length() < 3) continue;
	if (option[0]!='-' || option[1]!='-') continue;
	// OK, the option name starts with "--"; does it contain the '=' sign?
	idx = option.find("=", 2);
	if (idx == string::npos)
	{   // there is no '=' sign; we assume that this is a logical
	    // variable with the value 'true'
	    parameter_name = option.substr(2);
	    n = parameter_name.length();
	    for (j=0; j<n; j++)
		if (parameter_name[j]=='-') parameter_name[j] = ' ';
	    transform(parameter_name.begin(), parameter_name.end(),
		      parameter_name.begin(), (int(*)(int)) toupper);
	    p.value = "true";
	}
	else
	{
	    // convert the option name to a parameter name
	    parameter_name = option.substr(2,idx-2);
	    n = parameter_name.length();
	    for (j=0; j<n; j++)
		if (parameter_name[j]=='-') parameter_name[j] = ' ';
	    transform(parameter_name.begin(), parameter_name.end(),
		      parameter_name.begin(), (int(*)(int)) toupper);
	    // define the parameter value
	    p.value = option.substr(idx+1);
	    if (p.value.length()==0) continue;
	}
	// do we know a parameter with such a name?
	it = parameter_map.find(parameter_name);
	if (it != parameter_map.end())
	{ // just update the value
	    (it->second).value = p.value;
	}
	else
	{ // issue a warning and add the parameter to the map
	    cerr << "WARNING: the parameter '" << parameter_name
		 << "' is not known." << endl;
	    p.description = "";
	    parameter_map.insert(make_pair(parameter_name, p));
	}
    }
}

//-------------------------------------------------------------------------


Parameter& Parameters::operator() (const char* name)
{
    string parameter_name(name);
    map<string,Parameter>::iterator it;
    transform(parameter_name.begin(), parameter_name.end(),
	      parameter_name.begin(), (int(*)(int)) toupper);
    it = parameter_map.find(parameter_name);
    if (it == parameter_map.end())
    {
	string s = "The parameter '" + parameter_name +
	    "' has not been registered.";
	throw s;
    }
    return it->second;
}

//-------------------------------------------------------------------------


string Parameters::GetParameterValue(const char* name)
{
    string parameter_name(name);
    map<string,Parameter>::iterator it;
    transform(parameter_name.begin(), parameter_name.end(),
	      parameter_name.begin(), (int(*)(int)) toupper);
    it = parameter_map.find(parameter_name);
    if (it == parameter_map.end())
    {
	string s = "The parameter '" + parameter_name +
	    "' has not been registered.";
	throw s;
    }
    return (it->second).value;
}

//-------------------------------------------------------------------------


void Parameters::SetParameterValue(const char* name, string value)
{
    string parameter_name(name);
    map<string,Parameter>::iterator it;
    transform(parameter_name.begin(), parameter_name.end(),
	      parameter_name.begin(), (int(*)(int)) toupper);
    it = parameter_map.find(parameter_name);
    if (it == parameter_map.end())
    {
	string s = "The parameter '" + parameter_name +
	    "' has not been registered.";
	throw s;
    }
    (it->second).value = value;
}

//-------------------------------------------------------------------------


void Parameters::WriteParameterFile(ostream& os)
{
    map<string,Parameter>::iterator it;
    for (it=parameter_map.begin(); it!=parameter_map.end(); ++it)
    {
	PrintFormatted(os, (it->second).description, 75, "# ");
	os << (it->first) << " = " << (it->second).value << "\n\n";
    }
}

//-------------------------------------------------------------------------


void Parameters::PrintHelp()
{
    map<string,Parameter>::iterator it;
    string s;
    int i, n;

    s = "Each design is an ASCII file with two \
or three columns: the material name in the first column, the physical \
layer thickness in nanometres in the second one, and the optional sign '*' in \
the third column can instruct the code that the thickness of the layer may \
not be changed by optimisation. Before the code starts optimising or \
analysing the designs, it searches for the file 'parameters.txt' in the \
current directory and reads the parameters from this file. Each of the \
parameters in this file can be overridden by a command line option. For \
example, the option 'ANGLE OF INCIDENCE = 20' in the parameter file can \
be overridden by the command line option '--angle-of-incidence=0', and \
'VERBOSE = no' can be overridden by '--verbose=yes' or just \
'--verbose'. The valid command-line options are the following:";
    PrintFormatted(cout, s, 75, "");
    for (it=parameter_map.begin(); it!=parameter_map.end(); ++it)
    {
	s = "--" + it->first;
	transform(s.begin(), s.end(), s.begin(), (int(*)(int)) tolower);
	n = s.length();
	for (i=0; i<n; i++) if (s[i]==' ' || s[i]=='\t') s[i] = '-';
	cout << s << "\n";
    }
    exit(0);
}

//=========================================================================

Variables::Variables(Parameters& parameters) : parameters(parameters)
{
    target_dispersion = NULL;
}

Variables::~Variables()
{
    if (target_dispersion) delete target_dispersion;
}

AbstractDispersion& Variables::TargetDispersion()
{
    if (! target_dispersion)
    {
	int L;
	string str;
	Real omega0, lambda0, x;
	vector<Real> coefficients;
	str << parameters("TARGET DISPERSION");
	L = str.length();
	if (L==0)
	{
	    lambda0 << parameters("CENTRAL WAVELENGTH");
	    omega0 = 2.0*M_PI*SPEED_OF_LIGHT / lambda0;
	    coefficients.push_back(0);
	    TaylorDispersion* dispersion = new TaylorDispersion;
	    target_dispersion = dispersion;
	    dispersion->Set(omega0, coefficients);
	}
	else if ((str[0]=='"' && str[L-1]=='"') ||
	    (str[0]=='\'' && str[L-1]=='\''))
	{   // read the target dispersion from the file
	    str = str.substr(1, L-2);
	    FlexibleDispersion* dispersion = new FlexibleDispersion;
	    target_dispersion = dispersion;
	    dispersion->Set(str.c_str());
	}
	else
	{
	    lambda0 << parameters("CENTRAL WAVELENGTH");
	    omega0 = 2.0*M_PI*SPEED_OF_LIGHT / lambda0;
	    istringstream is(str);
	    while (is >> x) coefficients.push_back(x);
	    if (coefficients.size()==0)
	    {
		cerr << "ERROR: can't interpret '" << str
		     << "' as a list of whitespace-separated numbers;\n\
if you are trying to specify the file on the command line, try putting\n\
backslashes in front of the quotes" << endl;
		exit(1);
	    }
	    TaylorDispersion* dispersion = new TaylorDispersion;
	    target_dispersion = dispersion;
	    dispersion->Set(omega0, coefficients);
	}
    }
    return *target_dispersion;
}

void Variables::InitialiseProbePulse(const string& probe_pulse)
{
    int i, m, n;
    string str, s;
    vector<Real> wavelengths;
    vector<Real> intensity;
    vector<Real> phase;
    const int N=100; // number of points for the pulse, if it's not
		     // specified by a file
    int supergaussian_order, order;
    Real FWHM, spectral_width,lambda0, lambda_min, lambda_max;
    Real omega_min, omega_max, x, z, omega, lambda;
    istringstream is;

    n = probe_pulse.size();
    if (n==0) throw "ERROR in Variables::InitialiseProbePulse: \
the probe pulse is not defined";
    if ((probe_pulse[0]=='"' && probe_pulse[n-1]=='"') ||
	(probe_pulse[0]=='\'' && probe_pulse[n-1]=='\''))
    {  // read the spectrum from the file
	vector< vector<Real> > data;
	s = probe_pulse.substr(1,n-2);
	ifstream file(s.c_str());
	if (!file)
	{
	    str = "FAILED to open '" + s + "'";
	    throw str;
	}
	ReadTable(file, data);
	file.close();
	n = data.size();
	if (n < 1)
	{
	    cerr << "ERROR: no data have been read from '" << s.c_str()
		 << "'" << endl;
	    exit(1);
	}
	for (i=0; i<n; i++)
	{
	    m = data[i].size();
	    if (m < 2)
	    {
		str = "ERROR: can't parse the file '" + s + "'";
		throw str;
	    }
	    wavelengths.push_back(data[i][0]);
	    if (data[i][1] >= Real(0)) intensity.push_back(data[i][1]);
	    else intensity.push_back(0);
	    if (m>2) phase.push_back(data[i][2]);
	    else phase.push_back(0);
	}
    }
    else
    {
	// if the pulse is not specified by a file, it must be either
	// Gaussian or super-Gaussian
	is.str(probe_pulse);
	is >> str;
	transform(str.begin(), str.end(), str.begin(), (int(*)(int)) toupper);
	if (str=="GAUSSIAN") supergaussian_order = 1;
	else if (str=="SUPERGAUSSIAN")
	{
	    is >> supergaussian_order;
	    if (! is)
	    {
		str = "ERROR: can't obtain the order of the super-Gaussian \
pulse from\n\
PROBE PULSE = " + probe_pulse;
		throw str;
	    }
	}
	FWHM << parameters("PROBE PULSE FWHM");
	spectral_width = FWHM /
	    (2.0*pow(log(2.0), 0.5/supergaussian_order));
	lambda0 << parameters("PULSE CENTRAL WAVELENGTH");
	lambda_min << parameters("MINIMAL WAVELENGTH");
	lambda_max << parameters("MAXIMAL WAVELENGTH");
	omega_min = 2.0*M_PI*SPEED_OF_LIGHT / lambda_max;
	omega_max = 2.0*M_PI*SPEED_OF_LIGHT / lambda_min;
	for (i=0; i<N; i++)
	{
	    omega = omega_max - (omega_max - omega_min) * i/Real(N-1);
	    lambda = 2.0*M_PI*SPEED_OF_LIGHT / omega;
	    x = lambda - lambda0;
	    x = square(x / spectral_width);
	    z = 1.0;
	    for (order=0; order<supergaussian_order; order++) z *= x;
	    wavelengths.push_back(lambda);
	    intensity.push_back(exp(-z));
	    phase.push_back(0);
	}
    } // if ((probe_pulse[0]=='"' && probe_pulse[n-1]=='"') ||...
    probe_pulse_spectrum.TakeData(wavelengths, intensity);
    probe_pulse_phase.TakeData(wavelengths, phase);
}


VariableParameter& Variables::ProbePulseSpectrum()
{
    if (! probe_pulse_spectrum.isInitialised())
    {
	string str;
	str << parameters("PROBE PULSE");
	InitialiseProbePulse(str);
    }
    return probe_pulse_spectrum;
}

VariableParameter& Variables::ProbePulseSpectralPhase()
{
    if (! probe_pulse_phase.isInitialised())
    {
	string str;
	str << parameters("PROBE PULSE");
	InitialiseProbePulse(str);
    }
    return probe_pulse_phase;
}

VariableParameter& Variables::TargetReflectance()
{
    if (! target_reflectance.isInitialised())
    {
	string s;
	s << parameters("TARGET REFLECTANCE");
	if (target_reflectance.TakeData(s))
	{
	    cerr << "ERROR: can't get TARGET REFLECTANCE" << endl;
	    exit(1);
	}
    }
    return target_reflectance;
}


VariableParameter& Variables::ReflectanceTolerance()
{
    if (! reflectance_tolerance.isInitialised())
    {
	string s;
	s << parameters("REFLECTANCE TOLERANCE");
	if (reflectance_tolerance.TakeData(s))
	{
	    cerr << "ERROR: can't get REFLECTANCE TOLERANCE" << endl;
	    exit(1);
	}
    }
    return reflectance_tolerance;
}

VariableParameter& Variables::GDTolerance()
{
    if (! gd_tolerance.isInitialised())
    {
	string s;
	s << parameters("GD TOLERANCE");
	if (gd_tolerance.TakeData(s))
	{
	    cerr << "ERROR: can't get GD TOLERANCE" << endl;
	    exit(1);
	}
    }
    return gd_tolerance;
}

VariableParameter& Variables::GDDTolerance()
{
    if (! gdd_tolerance.isInitialised())
    {
	string s;
	s << parameters("GDD TOLERANCE");
	if (gdd_tolerance.TakeData(s))
	{
	    cerr << "ERROR: can't get GDD TOLERANCE" << endl;
	    exit(1);
	}
    }
    return gdd_tolerance;
}

VariableParameter& Variables::GDDReserve()
{
    if (! gdd_reserve.isInitialised())
    {
	string s;
	s << parameters("GDD RESERVE");
	if (gdd_reserve.TakeData(s))
	{
	    cerr << "ERROR: can't get GDD RESERVE" << endl;
	    exit(1);
	}
    }
    return gdd_reserve;
}

VariableParameter& Variables::GDReserve()
{
    if (! gd_reserve.isInitialised())
    {
	string s;
	s << parameters("GD RESERVE");
	if (gd_reserve.TakeData(s))
	{
	    cerr << "ERROR: can't get GD RESERVE\n" << endl;
	    exit(1);
	}
    }
    return gd_reserve;
}

VariableParameter& Variables::ReflectanceReserve()
{
    if (! reflectance_reserve.isInitialised())
    {
	string s;
	s << parameters("REFLECTANCE RESERVE");
	if (reflectance_reserve.TakeData(s))
	{
	    cerr << "ERROR: can't get REFLECTANCE RESERVE\n" << endl;
	    exit(1);
	}
    }
    return reflectance_reserve;
}

