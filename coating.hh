// The heart and soul of the project. Objects of the class 'Coating' can
// calculate the optical properties of coatings very very fast.

#ifndef __COATING_HH
#define __COATING_HH

#include "common.hh"
#include "material.hh"
#include "design.hh"
#include <vector>
#include <map>
#include <string>

class Layer;

//=========================================================================

class Coating
{

    // This class was designed for the optimisation algorithms, which
    // often change just one or a few layers. The populational algorithms
    // of global optimisation and the local-optimisation algorithms, which
    // use the gradient belong to this category. When all the layers change
    // their thicknesses, the evaluation of the coating properties takes
    // more time than could take, because we have to update the L- and
    // R- matrices.
    //
    // INSTRUCTIONS how to use the class:
    // ----------------------------------
    // 1) call "SetParameters" first;
    // 2) call "ImportDesign" when the parameters are initialised;
    // 3) after the first two steps you can call any member functions
    //    in any order.

public:
    Coating();
    Coating(const Coating& coating);
    const Coating& operator=(const Coating& coating);
    virtual ~Coating();

    typedef struct Parameters_tag
    {
	Tpolarisation polarisation;
	Real angle_of_incidence; // in radians
	string incidence_medium_name;
	string exit_medium_name;
	int N; // the number of frequencies to simulate
	Real omega_min, omega_max; // 1/fs
        /*
	  "omega_min" and "omega_max" are circular frequencies, which must be
	  given in '1/fs'; in order allow fast dispersion evaluation,
	  the coating properties are evaluated at "N" equidistant
	  frequencies between "omega_min" and "omega_max"
	*/
    } Parameters;

    void SetParameters(const Parameters& parameters,
		       MaterialRepository& material_repository);
    void GetParameters(Parameters& parameters) const;

    void ImportDesign(const Design& design,
		      MaterialRepository& material_repository);
    void ExportDesign(Design& design) const;

    int NumberOfLayers() const;
    void GetLayerThicknesses(Real* thicknesses) const;
    void SetLayerThicknesses(const Real* thicknesses);
    Real StackThickness() const; // the sum of all the layer thicknesses

    /* Some of the layer thicknesses can be fixed, so that they are not
     optimised. Each coating knows how many layers are subject to
     optimisation, it can return an array of thicknesses of
     these layers, and it can change layer thicknesses given
     such an array. */
    int NumberOfVariableLayers() const;
    void GetVariableThicknesses(Real* thicknesses) const;
    void SetVariableThicknesses(const Real* thicknesses);

    // the properties of the coating are calculated at certain frequencies
    int NumberOfFrequencies() const;
    Real FrequencyStep() const;
    void Frequencies(Real* f) const;
    void Wavelengths(Real* w) const;

    void Reflectivity(Complex* reflectivity);
    // computes the Fresnel coefficient of reflectivity at the registered
    // frequencies; the result is written to the array 'reflectivity'

    void TryLayerThickness(int n, Real thickness, Complex* reflectivity);
    /* Calculates reflectivity as if the variable layer thickness number 'n'
       were modified, without actually changing any layer thicknesses.  This
       function is very useful for gradient calculations, when we have to
       calculate the reflectivity for many designs, which differ from the
       basic one only in the thickness of exactly one layer */

    void Transmittance(Real* transmittance);
    // calculates the transmittance at the registered frequencies;
    // the result is written to the array 'transmittance'

    void Crossover(int crossover_point, const Coating& coating);
    /* This function is designated to the populational optimisation of
       coatings. "crossover_point" is the number of variable layers,
       which have to be taken from the current object, and the rest of
       the layers is taken from the object 'coating' passed as a parameter.
       The crossover could be implemented using 'GetVariableThicknesses'
       and 'SetVariableThicknesses', but 'Crossover' can save some time,
       because it doesn't have to evaluate any transfer matrices. The
       number of layers and the number of fixed layers must be the same
       in both objects, the corresponding layers must be made of the
       same materials. */

    void FieldAmplitude(Real penetration_depth, Complex* E, Complex* H,
			Complex* eta=NULL);
    /* This function returns the tangential electric and magnetic fields
       for every frequency at the given penetration depth. The fields
       are normalised so that at the boundary between the substrate and
       the first layer the tangential component of the amplitude of the
       electric field is equal to 1. If 'eta' is not NULL, it is initialised
       with the impedance of the medium at the given depth: H_t=eta*E_t
       for the light propagating towards the substrate. */

    void FieldAmplitude(int layer, Real z, Complex* E, Complex* H,
			Complex* eta=NULL);
    /* The same as the previous function, but restricted to a single layer.
       The 'z' coordinate is the distance to the surface of the layer,
       which is closer to the substrate. The zeroth layer is the one
       closest to the surface of the coating. If 'eta' is not NULL, it is
       initialised with the impedance of the medium of the given layer:
       H_t=eta*E_t for the light propagating towards the substrate. */

    void EFieldIntensity(const int N, const int Nx, const Real dx, 
	    Real** Z, const bool normalise=true);
    /* This function calculates the electric field intensity distribution inside
     * the layerstack. 'N' is the number of frequencies, 'Nx' is the number of 
     * penetration steps into the layers, 'dx' is the size of each step given 
     * in [nm]. If 'normalise' is true, the electric field intensity will be 
     * normalised to the incident field. In most cases this is convenient, thus 
     * the default value is 'true'.
     */

private:
    int N; // the number of frequencies to simulate
    int number_of_materials;
    Real omega_min, omega_max, omega_step; // 1/fs
    Tpolarisation polarisation;
    Real angle_of_incidence; // in radians
    string incidence_medium_name;
    string exit_medium_name;
    Real* temporary_thicknesses;
    vector<Layer> layers;

    // eta = n * cos(theta) for the TE wave (s-polarisation),
    //	     n / cos(theta) for the TM wave (p-polarisation);
    // 'eta' is a function of frequency
    vector<Real*> eta_re; // for each material
    vector<Real*> eta_im; // for each material

    // phase_factor = 2*pi*n*cos(theta)/lambda; if one multiplies it 
    // with the thickness of the layer, it will give the phase shift
    vector<Real*> phase_factor_re; // for each material
    vector<Real*> phase_factor_im; // for each material

    /* we have to be able to quickly find the right element in the
       arrays 'eta_re', 'eta_im', 'phase_factor_re', and
       'phase_factor_im' */
    vector<int> material_index; // maps a layer index to the material
    map<string,int> material_map; // maps material names to the index
    int incidence_medium_index, exit_medium_index;

    /* the following arrays are calculated each time when one or more layers
     change their thickness */
    Complex* B; // amplitude of the tangential magnetic field
    Complex* C; // amplitude of the tangential electric field

    void ChangeNumberOfFrequencies(int new_N);
    // deallocates memory if necessary and allocates memory if new_N>0

    void Copy(const Coating& coating);
    /* copies the given coating assuming that all the memory questions
       have already been resolved */

    int RegisterMaterial(const string& material_name,
			 MaterialRepository& material_repository);
    /* allocates the memory if necessary, calculates 'eta_re', 'eta_im',
       'phase_factor_re', 'phase_factor_im' for current parameters, and
       returns the index of the material
    */

    void UpdateTransferMatrix(int i);
    /* refills the transfer matrix of the given layer */

    void RefreshLayers();
    /* recalculates all the matrices owned by the Layer objects */

    void Calculate_B_and_C();
    /* assuming that all the transfer matrices as well as the 'L' and 'R'
       matrices have already been calculated, this function fills the
       arrays 'B' and 'C' */
};

//=========================================================================

void Reflectance(int N, Complex* reflectivity, Real* reflectance);
/* Calculates the reflectance given the reflectivity (which
   is supposed to be obtained by Coating::Reflectivity).
   "N" is the length of the arrays. */

void PhaseShift(int N, Complex* reflectivity, Real* phase_shift);
/* Calculates the phase shift given the reflectivity (which
   is supposed to be obtained by Coating::Reflectivity).
   "N" is the length of the arrays. The phase shift is returned in
   RADIANS, with the PI jumps are possibly eliminated. */

void PhaseToGD(int N, Real omega_step, Real* phase_shift, Real* GD);
/* Differentiates numerically the given array of the phase shift
   and writes (- d phi / d omega) to "GD" (if the phase shift is given
   in radians and the "omega_step" is measured in 1/fs, then the
   group delay will be in fs); notice that the first and the last
   elements of "GD" are not reliable. "N" is the length of the arrays. */

void PhaseToGDD(int N, Real omega_step, Real* phase_shift, Real* GDD);
/* Differentiates numerically the given array of the phase shift
   and writes (- d^2 phi / d omega^2) to "GDD" (if the phase shift is
   given in radians and the "omega_step" is measured in 1/fs, then the
   group delay dispersion will be in fs^2); notice that the first and
   the last elements of "GDD" are not reliable. "N" is the length of
   the arrays. */


#endif
