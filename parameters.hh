
#ifndef __PARAMETERS_HH
#define __PARAMETERS_HH

#include "common.hh"
#include "dispersion.hh"
#include <vector>
#include <string>
#include <map>
#include <iosfwd>

void PrintFormatted(ostream& os, const string& str,
		    int L, const string prefix);

//=========================================================================

/* There are some kinds of parameters than can be represented by either
   a scalar or by a table of (x, y) pairs or, possibly, by some
   functional dependence. The following classes handles such data. */

class AbstractParameter
{
public:
    AbstractParameter();
    virtual ~AbstractParameter();
    virtual Real operator()(Real x) const = 0;
};


class ConstantParameter : public AbstractParameter
{
public:
    ConstantParameter(Real value);
    ~ConstantParameter();
    ConstantParameter(const ConstantParameter&);
    ConstantParameter& operator=(const ConstantParameter&);
    Real operator()(Real x) const;
private:
    Real value;
};


class VariableParameter : public AbstractParameter
{
public:
    VariableParameter();
    ~VariableParameter();
    VariableParameter(const VariableParameter&);
    VariableParameter& operator=(const VariableParameter&);

    bool isInitialised();

    void TakeData(vector<Real> X, vector<Real> Y);
    void TakeData(const vector< pair<Real,Real> >&);

    int TakeData(const string& str);
    /* The parameter must be is either a number in the ASCII
       representation or a name of a file taken into quotes
       (both single and double quotes are allowed). In the first
       case, the function f(x) will always return the given
       number, in the second case the file will be read and
       f(x) will be determined by the linear interpolation. If no
       errors occur, the function returns zero. */
    
    Real operator()(Real x) const;
private:
    int N;
    vector< pair<Real,Real> > data;
};

//=========================================================================

class FlexibleDispersion : public AbstractDispersion
{  /* objects of this class are supposed to be initialised with
      some data read from a file */
public:
    FlexibleDispersion();
    FlexibleDispersion(const vector<Real>& wavelengths, 
		       const vector<Real>& GDD);
    ~FlexibleDispersion();
    void Set(const vector<Real>& wavelengths, const vector<Real>& GDD);
    void Set(const char* file_name);
    Real Phase(Real omega) const;
    Real GD(Real omega) const;
    Real GDD(Real omega) const;
private:
    VariableParameter phase_parameter;
    VariableParameter GD_parameter;
    VariableParameter GDD_parameter;
};

//=========================================================================

class Parameter
{
public:
    Parameter() {}
    ~Parameter() {}
    Parameter(const Parameter& p);
    const Parameter& operator=(const Parameter& p);
    string value;
    string description;
};

int operator << (int& x, const Parameter& p);
float operator << (float& x, const Parameter& p);
double operator << (double& x, const Parameter& p);
bool operator << (bool& x, const Parameter& p);
Tpolarisation operator << (Tpolarisation& x, const Parameter& p);
const string& operator << (string& x, const Parameter& p);
ostream& operator<< (ostream& o, const Parameter& p);


class Parameters
{
public:
    Parameters() {}
    Parameters(map<string,Parameter>& parameter_map);
    ~Parameters() {}

    void AddParameter(const char* name, const char* value,
		      const char* description);
    // the 'value' is considered to be the default value for the parameter

    void ParseFile(istream& is); // updates the parameters reading the file

    void ParseCommandLine(int argc, const char** argv);
    // updates the parameters parsing the command line; if the option '--help'
    // of '-h' is given on the command line, this function calls
    // PrintHelp() and terminates the program

    Parameter& operator() (const char* name);

    string GetParameterValue(const char* name); // returns the value
    void SetParameterValue(const char* name, string value);

    void WriteParameterFile(ostream& os);
    // send a template of the parameter file to the given stream

    void PrintHelp();
    // prints the parameter names from the parameter map and descriptions
    // of the parameters

private:
    map<string,Parameter> parameter_map;
};

//=========================================================================

class Variables
{
public:
    Variables(Parameters& parameters);
    ~Variables();
    AbstractDispersion& TargetDispersion();
    VariableParameter& TargetReflectance();
    VariableParameter& ReflectanceTolerance();
    VariableParameter& ReflectanceReserve();
    VariableParameter& GDTolerance();
    VariableParameter& GDReserve();
    VariableParameter& GDDTolerance();
    VariableParameter& GDDReserve();

    VariableParameter& ProbePulseSpectrum();
    /* the probe spectrum has to be a function of wavelength; if one
       wants the spectrum as a function of frequency, the spectral
       intensity has to be divided by omega^2 */

    VariableParameter& ProbePulseSpectralPhase();

private:
    Parameters& parameters;

    AbstractDispersion* target_dispersion;
    VariableParameter target_reflectance;
    VariableParameter reflectance_tolerance;
    VariableParameter reflectance_reserve;
    VariableParameter gd_tolerance;
    VariableParameter gd_reserve;
    VariableParameter gdd_tolerance;
    VariableParameter gdd_reserve;
    VariableParameter probe_pulse_spectrum;
    VariableParameter probe_pulse_phase;

    void InitialiseProbePulse(const string& probe_pulse);
};


#endif
