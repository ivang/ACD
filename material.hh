#ifndef __MATERIAL_HH
#define __MATERIAL_HH

#include "common.hh"
#include "interpolation.hh"
#include <string>
#include <vector>
#include <map>

//=========================================================================

class Material
{ // it knows how to read itself from the given data file, and how to use
  // the interpolation to find the refractive index at the given wavelength
public:
    Material(vector<string>* paths_to_materials);
    /* the data for each material is first searched in the current directory;
       if it's not found there, it is searched under the 'paths_to_materials'
    */
    virtual ~Material();
    Material(const Material& material);
    const Material& operator=(const Material& material);

    void ReadItself(const string& name);
    // 'name' is the name of the material, not the file

    Complex RefractiveIndex(Real lambda);
    bool ComplexIndex(Real lambda_min, Real lambda_max) const;
    // returns 'true' if the material is absorbing (gain is negative
    // absorption) within the specified wavelength range

    bool isMyName(const string& name) const;
    const string& MyName() const; // copies its name to the given string

protected:
    string _name;
    vector<string>* paths_to_materials;
    vector<Real> wavelengths;
    vector<Complex> refractive_indices;
    int N; // the number of data points read from the file
    AbstractInterpolator* interpolator;
    bool bandwidth_warning; // true, if the warning about the bandwidth
                            // violation has not been issued yet
    bool amplifying_medium; // true if some of the refractive indexes have
			    // a negative imaginary part (needed to detect
			    // the false amplification caused by the spline
			    // interpolation)
    bool amplification_warning; // true if the warning about the false
				// amplification has not been issued
};

//=========================================================================

class MaterialRepository
{ // Maintains a list of currently used materials
public:
    MaterialRepository(const vector<string>& paths_to_materials);
    /* the data for each material is first searched in the current directory;
       if it's not found there, it is searched under 'path_to_materials' */
    ~MaterialRepository();
    Complex RefractiveIndex(const string& material_name, Real lambda);

    bool ComplexIndex(const string& material_name,
		      Real lambda_min, Real lambda_max);
    // returns 'true' if the material is absorbing (gain is negative
    // absorption) within the specified wavelength range

    bool ComplexIndex(Real lambda_min, Real lambda_max);
    // returns 'true' if one of the already registered materials is
    // absorbing (gain is negative absorption) within the specified
    // wavelength range

    const Material& GetMaterial(const string& material_name);
    int NumberOfMaterials(); // returns the number of materials
			     // in the repository
private:
    vector<string> paths_to_materials;
    vector<Material> materials;
    map<string,int> index;

    int InsertMaterial(const string& name);
    // returns the position in 'materials', where the material has
    // been inserted
};

//=========================================================================
#endif
