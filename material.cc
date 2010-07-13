#include "material.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>

Material::Material(vector<string>* paths_to_materials) :
    paths_to_materials(paths_to_materials)
{
    bandwidth_warning = true;
    amplifying_medium = false;
    amplification_warning = true;
    interpolator = NULL;
    N = 0;
}

//-----------------------------------------

Material::Material(const Material& material)
{
    N = material.N;
    _name = material._name;
    paths_to_materials = material.paths_to_materials;
    wavelengths = material.wavelengths;
    refractive_indices = material.refractive_indices;
    if (N > 2)
    {
	interpolator =
	    new SplineInterpolator(wavelengths, refractive_indices);
    }
    else
    {
	interpolator =
	    new LinearInterpolator(wavelengths, refractive_indices);
    }
    bandwidth_warning = material.bandwidth_warning;
    amplifying_medium = material.amplifying_medium;
    amplification_warning = material.amplification_warning;
}

//-----------------------------------------

const Material& Material::operator=(const Material& material)
{
    N = material.N;
    _name = material._name;
    paths_to_materials = material.paths_to_materials;
    wavelengths = material.wavelengths;
    refractive_indices = material.refractive_indices;
    if (interpolator) delete interpolator;
    if (N > 2)
    {
	interpolator =
	    new SplineInterpolator(wavelengths, refractive_indices);
    }
    else
    {
	interpolator =
	    new LinearInterpolator(wavelengths, refractive_indices);
    }
    bandwidth_warning = material.bandwidth_warning;
    amplifying_medium = material.amplifying_medium;
    amplification_warning = material.amplification_warning;
    return *this;
}

//-----------------------------------------

Material::~Material()
{
    delete interpolator;
}

//-----------------------------------------

void Material::ReadItself(const string& name)
{
    int i, n, line_counter;
    const int max_string_length = 256;
    char str[max_string_length];
    istringstream str_stream;
    string file_name;
    Real wavelength, N_re, N_im;

    if (_name.size())
    { // the material has already read itself
	return;
    }
    // try in the current directory first
    file_name = name + ".dat";
    ifstream data_file(file_name.c_str());
    if (! data_file)
    {   // check the paths
	n = paths_to_materials->size();
	for (i=0; i<n; i++)
	{
	    file_name = (*paths_to_materials)[i] + name + ".dat";
	    data_file.clear();
	    data_file.open(file_name.c_str());
	    if (data_file) break;
	}
    }
    // if the file is not found, throw an exception
    if (! data_file)
    {
	string error_message = "no data found for the material " + name;
	throw error_message;
    }
    // read the file line by line
    line_counter = 0;
    while (data_file.getline(str, max_string_length))
    {
	line_counter++;
	// if the first non-blank character is '#', it's a comment string
	for (i=0; i<max_string_length && str[i]!='\0' && str[i]!='\n' &&
		      (str[i]==' ' || str[i]=='\t'); i++) {}
	if (i<max_string_length && 
	    (str[i]=='\0' || str[i]=='\n' || str[i]=='#')) continue;
	// is it an empty string?
	if (i==max_string_length || str[i]=='\0') continue;
	// try to interpret the string
	str_stream.str(string(str));
	str_stream.clear();
	if (!(str_stream >> wavelength >> N_re))
	{ // can't read two real numbers from the string
	    ostringstream os;
	    os << "can't parse line " << line_counter << " of file "
	       << file_name;
// 	    os << "\n" << str << "\nwavelength = " << wavelength 
// 	       << "\nN_re = " << N_re;
	    throw os.str();
	}
	wavelengths.push_back(wavelength);
	if ((str_stream >> N_im))
	    refractive_indices.push_back(Complex(N_re, -N_im));
	else refractive_indices.push_back(N_re);
    }
    data_file.close();
    N = wavelengths.size();
    if (N > 2)
    {
	interpolator =
	    new SplineInterpolator(wavelengths, refractive_indices);
    }
    else
    {
	interpolator =
	    new LinearInterpolator(wavelengths, refractive_indices);
    }
    // check whether the medium is amplifying
    for (i=0; i<N; i++)
    {
	if (imag(refractive_indices[i])>Real(0)) amplifying_medium = true;
    }
    _name = name;
}
//-----------------------------------------

bool Material::isMyName(const string& name) const
{
    return (name == _name);
}

const string& Material::MyName() const
{
    return _name;
}

//-----------------------------------------

bool Material::ComplexIndex(Real lambda_min, Real lambda_max) const
{
    for (int i=0; i<N; i++)
    {
	if (wavelengths[i]>=lambda_min && wavelengths[i]<=lambda_max &&
	    imag(refractive_indices[i])!=Real(0)) return true;
    }
    return false;
}

//-----------------------------------------

Complex Material::RefractiveIndex(Real lambda)
{
    Complex n;
    if (N>1 && (lambda < wavelengths[0] || lambda > wavelengths[N-1]) && 
	bandwidth_warning)
    {
	cerr << "WARNING: the refractive index data for the material '"
	     << _name << "' are specified\n\
only within the region " << wavelengths[0] <<"-"<< wavelengths[N-1]
	     << " nm; lambda = " << lambda << " nm is outside this range."
	     << endl;
	bandwidth_warning = false;
    }
    n = interpolator->operator()(lambda);
    // due to the spline interpolation an absorbing material may become
    // an amplifying medium at some wavelengths, which lead to artifacts
    // like a reflectance larger that 100%; in order to avoid it we
    // explicitly set gain to zero and issue a warning
    if (! amplifying_medium && imag(n)>Real(0))
    {
	n = real(n);
	if (amplification_warning)
	{
	    cerr << "\
WARNING: interpolation artifacts for the medium '" << _name << "'\
are detected\n\
(lambda = " << lambda << " nm). Make sure that the data points for the\n\
imaginary part of the refractive index lie on a smooth curve." << endl;
	    amplification_warning = false;
	}
    }
    return n;
}

//=========================================================================

MaterialRepository::
MaterialRepository(const vector<string>& paths)
{
    int i, n, L;
    n = paths.size();
    for (i=0; i<n; i++)
	if (paths[i].size()) paths_to_materials.push_back(paths[i]);
    // check if the paths end with the directory separator
    n = paths_to_materials.size();
    for (i=0; i<n; i++)
    {
	L = paths_to_materials[i].size();
	// if (L==0) continue;
	if (paths_to_materials[i][L-1] != directory_separator)
	    paths_to_materials[i] += directory_separator;
    }
}

MaterialRepository::~MaterialRepository()
{
}

int MaterialRepository::NumberOfMaterials()
{
    return materials.size();
}

int MaterialRepository::InsertMaterial(const string& name)
{
    int i;
    Material m(&paths_to_materials);
    m.ReadItself(name);
    i = materials.size();
    materials.push_back(m);
    index.insert(make_pair(name, i));
    return i;
}

Complex MaterialRepository::
RefractiveIndex(const string& material_name, Real lambda)
{
    int i;
    map<string,int>::const_iterator it = index.find(material_name);
    if (it==index.end()) i = InsertMaterial(material_name);
    else i = it->second;
    return materials[i].RefractiveIndex(lambda);
}

bool MaterialRepository::ComplexIndex(const string& material_name,
				      Real lambda_min, Real lambda_max)
{
    int i;
    map<string,int>::const_iterator it = index.find(material_name);
    if (it==index.end()) i = InsertMaterial(material_name);
    else i = it->second;
    return materials[i].ComplexIndex(lambda_min, lambda_max);
}

bool MaterialRepository::ComplexIndex(Real lambda_min, Real lambda_max)
{
    for (unsigned i=0; i<materials.size(); i++)
    {
	if (materials[i].ComplexIndex(lambda_min, lambda_max))
	    return true;
    }
    return false;
}

const Material& MaterialRepository::GetMaterial(const string& material_name)
{
    int i;
    map<string,int>::const_iterator it = index.find(material_name);
    if (it==index.end()) i = InsertMaterial(material_name);
    else i = it->second;
    return materials[i];
}

//=========================================================================

// int main()
// {
//     /* test "Material" */
//     Material m("./");
//     try
//     {
// 	m.ReadItself("SIO2");
// 	cout << "name: " <<  m.MyName() << endl;
// 	cout << "isMyName(" << m.MyName() << ") = "
// 	     << m.isMyName(m.MyName()) << endl;
// 	cout << "isMyName(xxx) -> "
// 	     << m.isMyName("xxx") << endl;
// 	cout << "n(800.0) = " << m.RefractiveIndex(800.0) << endl;
// 	cout << "ComplexIndex(100, 10000) -> " << m.ComplexIndex(100, 10000)
// 	     << endl;
//     }
//     catch(const string& error_message)
//     {
// 	cerr << error_message << endl;
// 	return 1;
//     }
//     catch(const char* error_message)
//     {
// 	cerr << error_message << endl;
// 	return 1;
//     }
//     /* test "MaterialRepository" */
//     try
//     {
// 	MaterialRepository r("./");
// 	cout << r.RefractiveIndex("AIR", 800.0) << endl;
// 	cout << r.RefractiveIndex("FUSI", 900.0) << endl;
// 	m = r.GetMaterial("FUSI");
// 	cout << "name: " <<  m.MyName() << endl;
// 	cout << r.NumberOfMaterials() << endl;
//     }
//     catch(const string& error_message)
//     {
// 	cerr << error_message << endl;
// 	return 1;
//     }
//     catch(const char* error_message)
//     {
// 	cerr << error_message << endl;
// 	return 1;
//     }
//     return 0;
// }
