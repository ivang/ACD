#include "design.hh"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <stdio.h>

Design::Design()
{
    no_of_layers = 0;
}

//--------------------------------------

Design::~Design()
{
    Clear();
}

//--------------------------------------

void Design::ReadFromFile(const string& file_name)
{
    int i, line_counter;
    const int max_string_length = 256;
    char str[max_string_length];
    istringstream str_stream;
    const string whitespace(" \t");
    string material, s;
    Real thickness;
    bool fixed;
    string::size_type pos;

    if (no_of_layers) Clear();
    ifstream data_file(file_name.c_str());
    // if the file is not found, try the default path
    if (! data_file)
    {
	string error_message = "can't open the design file " + file_name;
	throw error_message;
    }
    // read the file line by line
    line_counter = 0;
    no_of_layers = 0;
    while (data_file.getline(str, max_string_length))
    {
	line_counter++;
	// if the first non-blank character is '#', it's a comment string
	for (i=0; i<max_string_length && str[i]!='\0' && str[i]!='\n' &&
		 (str[i]==' ' || str[i]=='\t'); i++) {}
	if (i>=max_string_length) continue;
	if (str[i]=='\0' || str[i]=='#' || str[i]=='\n') continue;
	// try to interpret the string
	str_stream.str(string(str));
	str_stream.clear();
	if (!(str_stream >> material >> thickness))
	{
	    ostringstream os;
	    os << "can't parse line " << line_counter << " of file "
	       << file_name;
	    throw os.str();
	}
	// is the thickness of the design fixed?
	if ((str_stream >> s))
	{
	    pos = s.find_first_not_of(whitespace);
	    fixed = (pos!=string::npos && s[pos]=='*');
	}
	else fixed = false;
	// if the material is different from the material of the previous
	// layer of if it is the same, but one of the layers has a fixed
	// thickness, add the layer; otherwise increase the thickness of
	// the previous layer
	if (no_of_layers>0 && material==material_names[no_of_layers-1] &&
	    !fixed && !fixed_thicknesses[no_of_layers-1])
	{
	    thicknesses[no_of_layers-1] += thickness;
	}
	else
	{
	    material_names.push_back(material);
	    thicknesses.push_back(thickness);
	    fixed_thicknesses.push_back(fixed);
	    no_of_layers++;
	}
    }
    data_file.close();
    no_of_layers = material_names.size(); // to be sure
    name = file_name;
}

//--------------------------------------

void Design::SaveToFile(string file_name)
{
    int i, n, L;
    string backup_name;

    if (no_of_layers <= 0)
	throw string("ERROR in Design::SaveToFile: no_of_layers <= 0");
    if (! file_name.size()) file_name = name;
    // create a backup copy
    backup_name = file_name + ".bak";
    rename(file_name.c_str(), backup_name.c_str());
    // open the file
    ofstream file(file_name.c_str());
    if (! file)
    {
	string error_message = "can't save design '" + file_name + "'";
	throw error_message;
    }
    // find the longest material name
    L = 0;
    for (i=0; i<no_of_layers; i++)
    {
	n = material_names[i].size();
	if (L < n) L = n;
    }
    // save the design
    for (i=0; i<no_of_layers; i++)
    {
	file << setw(L) << left << material_names[i] <<" "
	     << setprecision(9) << setw(14) << right << thicknesses[i];
	if (fixed_thicknesses[i]) file << " *";
	file << "\n";
    }
    file.close();
}

//--------------------------------------

bool Design::isLike(const Design& other_design)
{
    if (no_of_layers != other_design.no_of_layers) return false;
    for (int i=0; i<no_of_layers; i++)
    {
	if (material_names[i] != other_design.material_names[i] ||
	    fixed_thicknesses[i] != other_design.fixed_thicknesses[i])
	    return false;
	if (fixed_thicknesses[i] && other_design.fixed_thicknesses[i] &&
	    thicknesses[i] != other_design.thicknesses[i]) return false;
    }
    return true;
}

//--------------------------------------

void Design::Clear()
{
    material_names.clear();
    thicknesses.clear();
    fixed_thicknesses.clear();
//    name.clear();
    no_of_layers = 0;
}

//--------------------------------------

void Design::AddLayer(const string& material_name, Real thickness,
		      bool fixed_thickness)
{
    material_names.push_back(material_name);
    thicknesses.push_back(thickness);
    fixed_thicknesses.push_back(fixed_thickness);
    no_of_layers++;
}

//--------------------------------------

const string& Design::MaterialName(int layer_index) const
{
    if (layer_index<0 || layer_index >= no_of_layers)
	throw "ERROR in Design::MaterialName: illegal layer_index";
    return material_names[layer_index];
}

//--------------------------------------

Real Design::LayerThickness(int layer_index) const
{
    if (layer_index<0 || layer_index >= no_of_layers)
	throw "ERROR in Design::LayerThickness: illegal layer_index";
    return thicknesses[layer_index];
}

//--------------------------------------

bool Design::FixedThickness(int layer_index) const
{
    if (layer_index<0 || layer_index >= no_of_layers)
	throw "ERROR in Design::FixedThickness: illegal layer_index";
    return fixed_thicknesses[layer_index];
}

Real Design::StackThickness() const
{
    Real stack_thickness = 0;
    for (int i=0; i<no_of_layers; i++) stack_thickness += thicknesses[i];
    return stack_thickness;
}

//-------------------------------------------------------------------------

void RenameDesignsToCWD(vector<Design>& designs)
{
    int i, j, k, n;
    ostringstream os;
    string name;
    bool designs_renamed;
    string::size_type pos;

    n = designs.size();
    // extract the path from each design name
    for (i=0; i<n; i++)
    {
	name = designs[i].name;
	pos = name.find_last_of(directory_separator);
	if (pos != string::npos) designs[i].name = name.substr(pos+1);
    }
    // make sure that all the designs have different names
    do
    {
	designs_renamed = false;
	for (i=0; i<n-1; i++)
	{
	    k = 0;
	    for (j=i+1; j<n; j++)
		if (designs[j].name == designs[i].name) k++;
	    if (k>0)
	    {   // rename the designs
		name = designs[i].name;
		designs[i].name += ".1";
		k = 2;
		for (j=i+1; j<n; j++)
		{
		    if (designs[j].name == name)
		    {
			os.str("");
			os << "." << k++;
			designs[j].name += os.str();
		    }
		}
		designs_renamed = true;
	    } // if (k>0)
	} // for (i=0; i<n-1; i++)
    } while(designs_renamed);
}


//=========================================================================

// int main()
// {
//     Design design;
//     try
//     {
// 	design.ReadFromFile("design.txt");
// 	cout << "name: " << design.name << endl;
// 	cout << "number of layers: " << design.NumberOfLayers() << endl;
// 	design.AddLayer("funny_material", 3.14, true);
// 	design.name = "design.xxx";
// 	design.SaveToFile();
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
