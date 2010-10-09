#include "coating.hh"
#include <cmath>

#ifdef DEBUG
#include <fstream>
#include <iomanip>
#endif

inline Complex ArcSin(Complex x)
{
    const Complex iu(Real(0),Real(1));
    if (imag(x)==0 && real(x)<=Real(1)) return asin(real(x));
    return -iu * log(iu*x + sqrt(Real(1) - x*x));
}

//=========================================================================

/* What makes the class 'Coating' special is the ability to calculate the
   properties of a coating very quickly, when just one or a few
   layers change its thicknesses. In order to achieve this, we
   store some information for each layer: the complex transfer
   matrix, as well as the left and the right products of transfer
   matrices (the matrix L is the product of all the transfer matrices
   from the first layer to current one excluding the current layer;
   the matrix R is the product of all the transfer matrices
   from the current layer to the last one excluding the current layer). */
class Layer
{
    friend class Coating;
public:
    Layer();
    Layer(int N, const string& material_name, Real thickness,
	  bool fixed_thickness);
    // N is the number of the frequencies, at which coating properties
    // have to be evaluated
    ~Layer();
    Layer(const Layer&);
    void Reset(int new_N, const string& new_material_name,
	       Real new_thickness, bool new_fixed_thickness);
    const Layer& operator=(const Layer& layer);

    /* if the refractive index of the material is real (that is, the
       material is neither absorbing nor has it gain), then the transfer
       matrix has real numbers on the main diagonal and pure complex
       off-diagonal elements; the following flags allow the code to
       choose the fastest way to calculate the L and R matrices */
    bool complex_A, complex_L, complex_R;

    bool L_matrix_updated;
    /* Coating::TryLayerThickness needs both, L- and R- matrices, but
       Coating::Calculate_B_and_C needs only the R-matrices. We can
       save some time if we only update the L-matrices, when they are
       required */

private:

    void UpdateL(const Layer& left_layer);
    // multiplies the L-matrix of the given layer by the transfer matrix
   
    void UpdateR(const Layer& right_layer);
    // multiplies the transfer matrix by the R-matrix of the given layer

    int N;
    string material_name;
    Real d; // layer thickness
    bool fixed_thickness; // true if the thickness should not be optimised
    Real *a11_re, *a11_im;
    Real *a12_re, *a12_im;
    Real *a21_re, *a21_im;
    // there is no need to store a22 since it is equal to a11
    Real *L11_re, *L11_im;
    Real *L12_re, *L12_im;
    Real *L21_re, *L21_im;
    Real *L22_re, *L22_im;
    Real *R11_re, *R11_im;
    Real *R12_re, *R12_im;
    Real *R21_re, *R21_im;
    Real *R22_re, *R22_im;
    void AllocateMemory();
    void DeallocateMemory();
    void Copy(const Layer& layer);
};

//-------------------------------------------------------------------------

Layer::Layer() : N(0), d(Real(0)), fixed_thickness(false)
{
    a11_re = NULL;
    a12_re = NULL;
    a21_re = NULL;
    L11_re = NULL;
    L12_re = NULL;
    L21_re = NULL;
    L22_re = NULL;
    R11_re = NULL;
    R12_re = NULL;
    R21_re = NULL;
    R22_re = NULL;
    a11_im = NULL;
    a12_im = NULL;
    a21_im = NULL;
    L11_im = NULL;
    L12_im = NULL;
    L21_im = NULL;
    L22_im = NULL;
    R11_im = NULL;
    R12_im = NULL;
    R21_im = NULL;
    R22_im = NULL;
    complex_A = true;
    complex_L = true;
    complex_R = true;
    L_matrix_updated = false;
}

//-------------------------------------------------------------------------

Layer::Layer(int N, const string& material_name, Real thickness,
	     bool fixed_thickness) :
    N(N), material_name(material_name), d(thickness),
    fixed_thickness(fixed_thickness)
{
    AllocateMemory();
    for (int i=0; i<N; i++)
    {
	a11_re[i] = 0.0;
	a11_im[i] = 0.0;
	a12_re[i] = 0.0;
	a12_im[i] = 0.0;
	a21_re[i] = 0.0;
	a21_im[i] = 0.0;
	L11_re[i] = 0.0;
	L11_im[i] = 0.0;
	L12_re[i] = 0.0;
	L12_im[i] = 0.0;
	L21_re[i] = 0.0;
	L21_im[i] = 0.0;
	L22_re[i] = 0.0;
	L22_im[i] = 0.0;
	R11_re[i] = 0.0;
	R11_im[i] = 0.0;
	R12_re[i] = 0.0;
	R12_im[i] = 0.0;
	R21_re[i] = 0.0;
	R21_im[i] = 0.0;
	R22_re[i] = 0.0;
	R22_im[i] = 0.0;
    }
    complex_A = true;
    complex_L = true;
    complex_R = true;
    L_matrix_updated = false;
}

//-------------------------------------------------------------------------

Layer::~Layer()
{
    DeallocateMemory();
}

//-------------------------------------------------------------------------

Layer::Layer(const Layer& layer) :
    N(layer.N), material_name(layer.material_name), d(layer.d),
    fixed_thickness(layer.fixed_thickness)
{
    AllocateMemory();
    Copy(layer);
}

//-------------------------------------------------------------------------

void Layer::Reset(int new_N, const string& new_material_name,
	   Real new_thickness, bool new_fixed_thickness)
{
    if (N != new_N)
    {
	DeallocateMemory();
	N = new_N;
	AllocateMemory();
    }
    material_name = new_material_name;
    d = new_thickness;
    fixed_thickness = new_fixed_thickness;
}

//-------------------------------------------------------------------------

const Layer& Layer::operator=(const Layer& layer)
{
    Reset(layer.N, layer.material_name, layer.d, layer.fixed_thickness);
    Copy(layer);
    return *this;
}

//-------------------------------------------------------------------------

void Layer::AllocateMemory()
{
    a11_re = new Real[N];
    a11_im = new Real[N];
    a12_re = new Real[N];
    a12_im = new Real[N];
    a21_re = new Real[N];
    a21_im = new Real[N];
    L11_re = new Real[N];
    L11_im = new Real[N];
    L12_re = new Real[N];
    L12_im = new Real[N];
    L21_re = new Real[N];
    L21_im = new Real[N];
    L22_re = new Real[N];
    L22_im = new Real[N];
    R11_re = new Real[N];
    R11_im = new Real[N];
    R12_re = new Real[N];
    R12_im = new Real[N];
    R21_re = new Real[N];
    R21_im = new Real[N];
    R22_re = new Real[N];
    R22_im = new Real[N];
}

//-------------------------------------------------------------------------

void Layer::DeallocateMemory()
{
    if (N > 0)
    {
	delete[] a11_re;
	delete[] a11_im;
	delete[] a12_re;
	delete[] a12_im;
	delete[] a21_re;
	delete[] a21_im;
	delete[] L11_re;
	delete[] L11_im;
	delete[] L12_re;
	delete[] L12_im;
	delete[] L21_re;
	delete[] L21_im;
	delete[] L22_re;
	delete[] L22_im;
	delete[] R11_re;
	delete[] R11_im;
	delete[] R12_re;
	delete[] R12_im;
	delete[] R21_re;
	delete[] R21_im;
	delete[] R22_re;
	delete[] R22_im;
    }
}

//-------------------------------------------------------------------------

void Layer::Copy(const Layer& layer)
{
    if (N != layer.N) throw "ERROR in Layer::Copy: N != layer.N";
    complex_A = layer.complex_A;
    complex_L = layer.complex_L;
    complex_R = layer.complex_R;
    L_matrix_updated = layer.L_matrix_updated;
    for (int i=0; i<N; i++)
    {
	a11_re[i] = layer.a11_re[i];
	a11_im[i] = layer.a11_im[i];
	a12_re[i] = layer.a12_re[i];
	a12_im[i] = layer.a12_im[i];
	a21_re[i] = layer.a21_re[i];
	a21_im[i] = layer.a21_im[i];
	L11_re[i] = layer.L11_re[i];
	L11_im[i] = layer.L11_im[i];
	L12_re[i] = layer.L12_re[i];
	L12_im[i] = layer.L12_im[i];
	L21_re[i] = layer.L21_re[i];
	L21_im[i] = layer.L21_im[i];
	L22_re[i] = layer.L22_re[i];
	L22_im[i] = layer.L22_im[i];
	R11_re[i] = layer.R11_re[i];
	R11_im[i] = layer.R11_im[i];
	R12_re[i] = layer.R12_re[i];
	R12_im[i] = layer.R12_im[i];
	R21_re[i] = layer.R21_re[i];
	R21_im[i] = layer.R21_im[i];
	R22_re[i] = layer.R22_re[i];
	R22_im[i] = layer.R22_im[i];
    }
}

//-------------------------------------------------------------------------

void Layer::UpdateL(const Layer& l)
{
    int i;
    if (N != l.N) throw "ERROR in Layer::UpdateL: N != l.N";
    if (! l.L_matrix_updated) throw "ERROR in Layer::UpdateL: \
the L-matrix of the given layer is not updated";
    complex_L = (l.complex_A || l.complex_L);
    if (l.complex_A)
    {
	if (l.complex_L)
	{
	    for (i=0; i<N; i++)
	    {
		L11_re[i] = l.L11_re[i]*l.a11_re[i] -
		    l.L11_im[i]*l.a11_im[i] + l.L12_re[i]*l.a21_re[i] -
		    l.L12_im[i]*l.a21_im[i];
		L11_im[i] = l.L11_re[i]*l.a11_im[i] +
		    l.L11_im[i]*l.a11_re[i] + l.L12_re[i]*l.a21_im[i] +
		    l.L12_im[i]*l.a21_re[i];
		L12_re[i] = l.L11_re[i]*l.a12_re[i] -
		    l.L11_im[i]*l.a12_im[i] + l.L12_re[i]*l.a11_re[i] -
		    l.L12_im[i]*l.a11_im[i];
		L12_im[i] = l.L11_re[i]*l.a12_im[i] +
		    l.L11_im[i]*l.a12_re[i] + l.L12_re[i]*l.a11_im[i] +
		    l.L12_im[i]*l.a11_re[i];
		L21_re[i] = l.L21_re[i]*l.a11_re[i] -
		    l.L21_im[i]*l.a11_im[i] + l.L22_re[i]*l.a21_re[i] -
		    l.L22_im[i]*l.a21_im[i];
		L21_im[i] = l.L21_re[i]*l.a11_im[i] +
		    l.L21_im[i]*l.a11_re[i] + l.L22_re[i]*l.a21_im[i] +
		    l.L22_im[i]*l.a21_re[i];
		L22_re[i] = l.L21_re[i]*l.a12_re[i] -
		    l.L21_im[i]*l.a12_im[i] + l.L22_re[i]*l.a11_re[i] -
		    l.L22_im[i]*l.a11_im[i];
		L22_im[i] = l.L21_re[i]*l.a12_im[i] +
		    l.L21_im[i]*l.a12_re[i] + l.L22_re[i]*l.a11_im[i] +
		    l.L22_im[i]*l.a11_re[i];
	    }
	}
	else
	{
	    for (i=0; i<N; i++)
	    {
		L11_re[i] = l.L11_re[i]*l.a11_re[i] - l.L12_im[i]*l.a21_im[i];
		L11_im[i] = l.L11_re[i]*l.a11_im[i] + l.L12_im[i]*l.a21_re[i];
		L12_re[i] = l.L11_re[i]*l.a12_re[i] - l.L12_im[i]*l.a11_im[i];
		L12_im[i] = l.L11_re[i]*l.a12_im[i] + l.L12_im[i]*l.a11_re[i];
		L21_re[i] = l.L22_re[i]*l.a21_re[i] - l.L21_im[i]*l.a11_im[i];
		L21_im[i] = l.L21_im[i]*l.a11_re[i] + l.L22_re[i]*l.a21_im[i];
		L22_re[i] = l.L22_re[i]*l.a11_re[i] - l.L21_im[i]*l.a12_im[i];
		L22_im[i] = l.L21_im[i]*l.a12_re[i] + l.L22_re[i]*l.a11_im[i];
	    }
	} // if (l.complex_L)
    }
    else // l.complex_A==false
    {
	if (l.complex_L)
	{
	    for (i=0; i<N; i++)
	    {
		L11_re[i] = l.L11_re[i]*l.a11_re[i] - l.L12_im[i]*l.a21_im[i];
		L11_im[i] = l.L11_im[i]*l.a11_re[i] + l.L12_re[i]*l.a21_im[i];
		L12_re[i] = l.L12_re[i]*l.a11_re[i] - l.L11_im[i]*l.a12_im[i];
		L12_im[i] = l.L11_re[i]*l.a12_im[i] + l.L12_im[i]*l.a11_re[i];
		L21_re[i] = l.L21_re[i]*l.a11_re[i] - l.L22_im[i]*l.a21_im[i];
		L21_im[i] = l.L21_im[i]*l.a11_re[i] + l.L22_re[i]*l.a21_im[i];
		L22_re[i] = l.L22_re[i]*l.a11_re[i] - l.L21_im[i]*l.a12_im[i];
		L22_im[i] = l.L21_re[i]*l.a12_im[i] + l.L22_im[i]*l.a11_re[i];
	    }
	}
	else
	{
	    for (i=0; i<N; i++)
	    {
		L11_re[i] = l.L11_re[i]*l.a11_re[i] - l.L12_im[i]*l.a21_im[i];
		L11_im[i] = 0.0;
		L12_re[i] = 0.0;
		L12_im[i] = l.L11_re[i]*l.a12_im[i] + l.L12_im[i]*l.a11_re[i];
		L21_re[i] = 0.0;
		L21_im[i] = l.L21_im[i]*l.a11_re[i] + l.L22_re[i]*l.a21_im[i];
		L22_re[i] = l.L22_re[i]*l.a11_re[i] - l.L21_im[i]*l.a12_im[i];
		L22_im[i] = 0.0;
	    }
	} // if (l.complex_L)
    } // if (l.complex_A)
    L_matrix_updated = true;
}

//-------------------------------------------------------------------------

void Layer::UpdateR(const Layer& l)
{
    int i;
    if (N != l.N) throw "ERROR in Layer::UpdateR: N != l.N";
    complex_R = (l.complex_A || l.complex_R);
    if (l.complex_A)
    {
	if (l.complex_R)
	{
	    for (i=0; i<N; i++)
	    {
		R11_re[i] = l.a11_re[i]*l.R11_re[i] -
		    l.a11_im[i]*l.R11_im[i] + l.a12_re[i]*l.R21_re[i] -
		    l.a12_im[i]*l.R21_im[i];
		R11_im[i] = l.a11_re[i]*l.R11_im[i] +
		    l.a11_im[i]*l.R11_re[i] + l.a12_re[i]*l.R21_im[i] +
		    l.a12_im[i]*l.R21_re[i];
		R12_re[i] = l.a11_re[i]*l.R12_re[i] -
		    l.a11_im[i]*l.R12_im[i] + l.a12_re[i]*l.R22_re[i] -
		    l.a12_im[i]*l.R22_im[i];
		R12_im[i] = l.a11_re[i]*l.R12_im[i] +
		    l.a11_im[i]*l.R12_re[i] + l.a12_re[i]*l.R22_im[i] +
		    l.a12_im[i]*l.R22_re[i];
		R21_re[i] = l.a21_re[i]*l.R11_re[i] -
		    l.a21_im[i]*l.R11_im[i] + l.a11_re[i]*l.R21_re[i] -
		    l.a11_im[i]*l.R21_im[i];
		R21_im[i] = l.a21_re[i]*l.R11_im[i] +
		    l.a21_im[i]*l.R11_re[i] + l.a11_re[i]*l.R21_im[i] +
		    l.a11_im[i]*l.R21_re[i];
		R22_re[i] = l.a21_re[i]*l.R12_re[i] -
		    l.a21_im[i]*l.R12_im[i] + l.a11_re[i]*l.R22_re[i] -
		    l.a11_im[i]*l.R22_im[i];
		R22_im[i] = l.a21_re[i]*l.R12_im[i] +
		    l.a21_im[i]*l.R12_re[i] + l.a11_re[i]*l.R22_im[i] +
		    l.a11_im[i]*l.R22_re[i];
	    }
	}
	else
	{
	    for (i=0; i<N; i++)
	    {
		R11_re[i] = l.a11_re[i]*l.R11_re[i] - l.a12_im[i]*l.R21_im[i];
		R11_im[i] = l.a11_im[i]*l.R11_re[i] + l.a12_re[i]*l.R21_im[i];
		R12_re[i] = l.a12_re[i]*l.R22_re[i] - l.a11_im[i]*l.R12_im[i];
		R12_im[i] = l.a11_re[i]*l.R12_im[i] + l.a12_im[i]*l.R22_re[i];
		R21_re[i] = l.a21_re[i]*l.R11_re[i] - l.a11_im[i]*l.R21_im[i];
		R21_im[i] = l.a21_im[i]*l.R11_re[i] + l.a11_re[i]*l.R21_im[i];
		R22_re[i] = l.a11_re[i]*l.R22_re[i] - l.a21_im[i]*l.R12_im[i];
		R22_im[i] = l.a21_re[i]*l.R12_im[i] + l.a11_im[i]*l.R22_re[i];
	    }
	} // if (l.complex_R)
    }
    else // l.complex_A==false
    {
	if (l.complex_R)
	{
	    for (i=0; i<N; i++)
	    {
		R11_re[i] = l.a11_re[i]*l.R11_re[i] - l.a12_im[i]*l.R21_im[i];
		R11_im[i] = l.a11_re[i]*l.R11_im[i] + l.a12_im[i]*l.R21_re[i];
		R12_re[i] = l.a11_re[i]*l.R12_re[i] - l.a12_im[i]*l.R22_im[i];
		R12_im[i] = l.a11_re[i]*l.R12_im[i] + l.a12_im[i]*l.R22_re[i];
		R21_re[i] = l.a11_re[i]*l.R21_re[i] - l.a21_im[i]*l.R11_im[i];
		R21_im[i] = l.a21_im[i]*l.R11_re[i] + l.a11_re[i]*l.R21_im[i];
		R22_re[i] = l.a11_re[i]*l.R22_re[i] - l.a21_im[i]*l.R12_im[i];
		R22_im[i] = l.a21_im[i]*l.R12_re[i] + l.a11_re[i]*l.R22_im[i];
	    }
	}
	else
	{
	    for (i=0; i<N; i++)
	    {
		R11_re[i] = l.a11_re[i]*l.R11_re[i] - l.a12_im[i]*l.R21_im[i];
		R11_im[i] = 0.0;
		R12_re[i] = 0.0;
		R12_im[i] = l.a11_re[i]*l.R12_im[i] + l.a12_im[i]*l.R22_re[i];
		R21_re[i] = 0.0;
		R21_im[i] = l.a21_im[i]*l.R11_re[i] + l.a11_re[i]*l.R21_im[i];
		R22_re[i] = l.a11_re[i]*l.R22_re[i] - l.a21_im[i]*l.R12_im[i];
		R22_im[i] = 0.0;
	    }
	} // if (l.complex_R)
    } // if (l.complex_A)
}

//=========================================================================

Coating::Coating()
{
    N = 0;
    number_of_materials = 0;
    omega_min = 0.0;
    omega_max = 0.0;
    omega_step = 0.0;
    polarisation = TM;
    angle_of_incidence = 0.0;
    incidence_medium_index = 0;
    exit_medium_index = 0;
    B = NULL;
    C = NULL;
    E_0 = NULL;
    H_0 = NULL;
    eta_0 = NULL;
    temporary_thicknesses = NULL;
}

//-------------------------------------------------------------------------

void Coating::ChangeNumberOfFrequencies(int new_N)
{
    int i;
    if (N==new_N) return;
    if (new_N == 0 || N > 0)
    {
	for (i=0; i<number_of_materials; i++)
	{
	    delete[] eta_re[i];
	    delete[] eta_im[i];
	    delete[] phase_factor_re[i];
	    delete[] phase_factor_im[i];
	}
	delete[] B;
	delete[] C;
	delete[] E_0;
	delete[] H_0;
	delete[] eta_0;
	for (i=0; i<NumberOfLayers(); i++)
	{
	    layers[i].DeallocateMemory();
	    layers[i].N = 0;
	}
    }
    N = new_N;
    if (N > 0)
    {
	for (i=0; i<number_of_materials; i++)
	{
	    eta_re[i] = new Real[N];
	    eta_im[i] = new Real[N];
	    phase_factor_re[i] = new Real[N];
	    phase_factor_im[i] = new Real[N];
	}
	B = new Complex[N];
	C = new Complex[N];
	for (i=0; i<NumberOfLayers(); i++)
	{
	    layers[i].N = N;
	    layers[i].AllocateMemory();
	}
	E_0 = new Complex[N];
	H_0 = new Complex[N];
	eta_0 = new Complex[N];
    }
}

//-------------------------------------------------------------------------

void Coating::Copy(const Coating& coating)
{
    int n;
    // copy the arrays
    for (int i=0; i<number_of_materials; i++)
    {
	for (n=0; n<N; n++)
	{
	    eta_re[i][n] = coating.eta_re[i][n];
	    eta_im[i][n] = coating.eta_im[i][n];
	    phase_factor_re[i][n] = coating.phase_factor_re[i][n];
	    phase_factor_im[i][n] = coating.phase_factor_im[i][n];
	}
    }
    for (n=0; n<N; n++)
    {
	B[n] = coating.B[n];
	C[n] = coating.C[n];
    }

    // copy the variables and simple containers
    omega_min = coating.omega_min;
    omega_max = coating.omega_max;
    omega_step = coating.omega_step;
    polarisation = coating.polarisation;
    angle_of_incidence = coating.angle_of_incidence;
    incidence_medium_name = coating.incidence_medium_name;
    exit_medium_name = coating.exit_medium_name;
    layers = coating.layers;
    incidence_medium_index = coating.incidence_medium_index;
    exit_medium_index = coating.exit_medium_index;
    material_index = coating.material_index;
    material_map = coating.material_map;
	
    // calculate the incident electric field
    FieldAmplitude(-1, E_0, H_0, eta_0);
}

//-------------------------------------------------------------------------

Coating::Coating(const Coating& coating)
{
    // memory management
    number_of_materials = coating.number_of_materials;
    eta_re.resize(number_of_materials);
    eta_im.resize(number_of_materials);
    phase_factor_re.resize(number_of_materials);
    phase_factor_im.resize(number_of_materials);
    temporary_thicknesses = new Real[coating.NumberOfLayers()];
    N = 0;
    ChangeNumberOfFrequencies(coating.N);
    Copy(coating);
}

//-------------------------------------------------------------------------

const Coating& Coating::operator=(const Coating& coating)
{
    int n;
    // memory management
    n = coating.NumberOfLayers();
    if (n != NumberOfLayers()) delete[] temporary_thicknesses;
    temporary_thicknesses = new Real[n];
    if (number_of_materials != coating.number_of_materials)
    {
	n = coating.number_of_materials;
	ChangeNumberOfFrequencies(0);
	eta_re.resize(n);
	eta_im.resize(n);
	phase_factor_re.resize(n);
	phase_factor_im.resize(n);
	number_of_materials = coating.number_of_materials;
    }
    ChangeNumberOfFrequencies(coating.N);
    Copy(coating);
    return *this;
}

//-------------------------------------------------------------------------

Coating::~Coating()
{
    ChangeNumberOfFrequencies(0);
    delete[] temporary_thicknesses;
}

//-------------------------------------------------------------------------

void Coating::SetParameters(const Parameters& parameters,
			    MaterialRepository& material_repository)
{
    polarisation = parameters.polarisation;
    angle_of_incidence = parameters.angle_of_incidence;
    incidence_medium_name = parameters.incidence_medium_name;
    exit_medium_name = parameters.exit_medium_name;
    omega_min = parameters.omega_min;
    omega_max = parameters.omega_max;
    ChangeNumberOfFrequencies(parameters.N);
    // refill the arrays, which reflect the material properties
    map<string,int>::const_iterator it;
    for (it=material_map.begin(); it!=material_map.end(); it++)
    {
	RegisterMaterial(it->first, material_repository);
    }
    // check if we have registered the incidence and exit media
    it = material_map.find(incidence_medium_name);
    if (it==material_map.end())
    {
	incidence_medium_index =
	    RegisterMaterial(incidence_medium_name, material_repository);
    }
    else incidence_medium_index = it->second;
    it = material_map.find(exit_medium_name);
    if (it==material_map.end())
    {
	exit_medium_index =
	    RegisterMaterial(exit_medium_name, material_repository);
    }
    else exit_medium_index = it->second;
    if (NumberOfLayers() > 0)
    {
	// refill the arrays of the Layer objects
	RefreshLayers();
	// calculate the tangential magnetic and electric fields
	Calculate_B_and_C();
	// calculate the incident electric field
	//FieldAmplitude(-1, E_0, H_0, eta_0);
    }
}

//-------------------------------------------------------------------------

void Coating::GetParameters(Parameters& parameters) const
{
    parameters.polarisation = polarisation;
    parameters.angle_of_incidence = angle_of_incidence;
    parameters.incidence_medium_name = incidence_medium_name;
    parameters.exit_medium_name = exit_medium_name;
    parameters.N = N;
    parameters.omega_min = omega_min;
    parameters.omega_max = omega_max;
}

//-------------------------------------------------------------------------

void Coating::ImportDesign(const Design& design,
			   MaterialRepository& material_repository)
{
    int i, j, n;
    Real lambda_min, lambda_max;
    map<string,int>::const_iterator it;

    if (omega_min==Real(0) || omega_max==Real(0))
	throw "ERROR in Coating::ImportDesign";
    lambda_min = 2.0*M_PI*SPEED_OF_LIGHT / omega_max;
    lambda_max = 2.0*M_PI*SPEED_OF_LIGHT / omega_min;
    n = design.NumberOfLayers();
    layers.resize(n);
    material_index.resize(n);
    if (temporary_thicknesses) delete[] temporary_thicknesses;
    temporary_thicknesses = new Real[n];
    /* For historical reasons the layers are inserted in the reverse order.
       This is necessary, because the transfer matrix is equal to the
       matrix product M1*M2*...*Mn, where M1 corresponds to the last layer,
       and Mn corresponds to the first layer in the layer stack, where the
       layers are counted from the substrate. */
    for (i=0, j=n-1; i<n; i++, j--)
    {
	// maybe the layer is made from a material we don't know yet
	it = material_map.find(design.MaterialName(i));
	if (it==material_map.end())
	{
	    material_index[j] =
		RegisterMaterial(design.MaterialName(i), material_repository);
	}
	else // we know the material, but we have to update the material_index
	{
	    material_index[j] = it->second;
	}
	// change the parameters of the layer
	layers[j].Reset(N, design.MaterialName(i),
			abs(design.LayerThickness(i)),
			design.FixedThickness(i));
    }
    // refill the matrices layers know about
    RefreshLayers();
    // calculate the tangential magnetic and electric fields
    Calculate_B_and_C();
    // calculate the incident electric field
    FieldAmplitude(-1, E_0, H_0, eta_0);
}

//-------------------------------------------------------------------------

void Coating::ExportDesign(Design& design) const
{
    int i, n;
    design.Clear();
    n = NumberOfLayers();
    for (i=0; i<n; i++)
    {
	design.AddLayer(layers[n-i-1].material_name, layers[n-i-1].d,
			layers[n-i-1].fixed_thickness);
    }
}

//-------------------------------------------------------------------------

int Coating::NumberOfLayers() const
{
    return layers.size();
}

//-------------------------------------------------------------------------

void Coating::GetLayerThicknesses(Real* thicknesses) const
{
    int i, n;
    n = NumberOfLayers();
    for (i=0; i<n; i++) thicknesses[i] = layers[i].d;
}

//-------------------------------------------------------------------------

void Coating::SetLayerThicknesses(const Real* thicknesses)
{
    int i, n;
    Real original_thickness;
    n = NumberOfLayers();
    if (n==0) throw "ERROR in Coating::SetLayerThicknesses";
    // recalculate the transfer matrices, but leave the layer thicknesses
    // as they are
    for (i=0; i<n; i++)
    {
	original_thickness = layers[i].d;
	if (thicknesses[i] != original_thickness)
	{
	    layers[i].d = abs(thicknesses[i]);
	    UpdateTransferMatrix(i);
	    layers[i].d = original_thickness;
	}
    }
    // mark the L matrices, which are out of date
    for (i=0; i<n; i++) if (layers[i].d != abs(thicknesses[i])) break;
    for (i++; i<n; i++) layers[i].L_matrix_updated = false;
//    for (i++; i<n; i++) layers[i].UpdateL(layers[i-1]);
    // update the R matrices
    for (i=n-1; i>=0; i--) if (layers[i].d != abs(thicknesses[i])) break;
    for (i--; i>=0; i--) layers[i].UpdateR(layers[i+1]);
    // save the new layer thicknesses
    for (i=0; i<n; i++) layers[i].d = abs(thicknesses[i]);
    // calculate the tangential magnetic and electric fields
    Calculate_B_and_C();
    // calculate the incident electric field
    FieldAmplitude(-1, E_0, H_0, eta_0);
}

//-------------------------------------------------------------------------

Real Coating::StackThickness() const
{
    int i, n;
    Real stack_thickness = 0;
    n = NumberOfLayers();
    for (i=0; i<n; i++) stack_thickness += layers[i].d;
    return stack_thickness;
}

//-------------------------------------------------------------------------

int Coating::NumberOfVariableLayers() const
{
    int i, j, n;
    n = NumberOfLayers();
    for (i=0, j=0; i<n; i++)
    {
	if (! layers[i].fixed_thickness) j++;
    }
    return j;
}

//-------------------------------------------------------------------------

void Coating::GetVariableThicknesses(Real* thicknesses) const
{
    int i, j, n;
    n = NumberOfLayers();
    for (i=0, j=0; i<n; i++)
    {
	if (! layers[i].fixed_thickness) thicknesses[j++] = layers[i].d;
    }
}

//-------------------------------------------------------------------------

void Coating::SetVariableThicknesses(const Real* thicknesses)
{
    int i, j, n;
    n = NumberOfLayers();
    if (n==0) throw "ERROR in Coating::SetVariableThicknesses";
    for (i=0, j=0; i<n; i++)
    {
	if (layers[i].fixed_thickness)
	    temporary_thicknesses[i] = layers[i].d;
	else temporary_thicknesses[i] = thicknesses[j++];
    }
    SetLayerThicknesses(temporary_thicknesses);
}

//-------------------------------------------------------------------------

int Coating::NumberOfFrequencies() const
{
    return N;
}

//-------------------------------------------------------------------------

Real Coating::FrequencyStep() const
{
    if (N > 0) return (omega_max - omega_min) / Real(N-1);
    return 0.0;
}

//-------------------------------------------------------------------------

void Coating::Frequencies(Real* f) const
{
    for (int i=0; i<N; i++)
	f[i] = omega_min + (omega_max - omega_min) * i/Real(N-1);
}

//-------------------------------------------------------------------------

void Coating::Wavelengths(Real* w) const
{
    Real omega;
    for (int i=0; i<N; i++)
    {
	omega = omega_min + (omega_max - omega_min) * i/Real(N-1);
	w[i] = 2.0*M_PI*SPEED_OF_LIGHT / omega;
    }
}

//-------------------------------------------------------------------------

void Coating::Reflectivity(Complex* reflectivity)
{
    Complex Y;
    Complex eta0;
    for (int i=0; i<N; i++)
    {
	Y = C[i] / B[i]; // the optical admittance
	eta0 = Complex(eta_re[incidence_medium_index][i],
		       eta_im[incidence_medium_index][i]);
	reflectivity[i] = (eta0 - Y) / (eta0 + Y);
    }
}

//-------------------------------------------------------------------------

void Coating::Transmittance(Real* transmittance)
{
    Complex eta0, eta1;
    for (int i=0; i<N; i++)
    {
	eta0 = Complex(eta_re[incidence_medium_index][i],
		       eta_im[incidence_medium_index][i]);
	eta1 = Complex(eta_re[exit_medium_index][i],
		       eta_im[exit_medium_index][i]);
	transmittance[i] = 4.0*real(eta0)*real(eta1) / norm(eta0*B[i] + C[i]);
    }
}

//-------------------------------------------------------------------------

void Coating::TryLayerThickness(int n, Real thickness, Complex* reflectivity)
{
    int i, j, number_of_layers;
    Real original_thickness;
    Complex Y; //, eta;
    Real A11_re, A11_im, A12_re, A12_im, A21_re, A21_im;
    Real L11_re, L11_im, L12_re, L12_im, L21_re, L21_im, L22_re, L22_im;
    Real R11_re, R11_im, R12_re, R12_im, R21_re, R21_im, R22_re, R22_im;
    Real M11_re, M11_im, M12_re, M12_im, M21_re, M21_im, M22_re, M22_im;
    Real eta_real, eta_imag;
    Complex eta0;
    bool complex_A;

    thickness = abs(thickness);

    // identify the layer, which is supposed to change its thickness
    number_of_layers = NumberOfLayers();
    if (number_of_layers==0) throw "ERROR in Coating::TryLayerThickness";
    for (i=0, j=0; i<number_of_layers && j<n; i++)
    {
	if (! layers[i].fixed_thickness) j++;
    }
    // it's the i-th layer
    original_thickness = layers[i].d;
    layers[i].d = thickness;
    UpdateTransferMatrix(i);
    // update the L-matrices if necessary
    for (j=1; j<=i; j++) if (! layers[j].L_matrix_updated) break;
    for (; j<=i; j++) layers[j].UpdateL(layers[j-1]);
    // calculate the reflectivity
    complex_A = layers[i].complex_A;
    for (j=0; j<N; j++)
    {
	// prepare the complex transfer-, L-, and R-matrices of the layer
	A11_re = layers[i].a11_re[j];
	A11_im = layers[i].a11_im[j];
	A12_re = layers[i].a12_re[j];
	A12_im = layers[i].a12_im[j];
	A21_re = layers[i].a21_re[j];
	A21_im = layers[i].a21_im[j];
	L11_re = layers[i].L11_re[j];
	L11_im = layers[i].L11_im[j];
	L12_re = layers[i].L12_re[j];
	L12_im = layers[i].L12_im[j];
	L21_re = layers[i].L21_re[j];
	L21_im = layers[i].L21_im[j];
	L22_re = layers[i].L22_re[j];
	L22_im = layers[i].L22_im[j];
	R11_re = layers[i].R11_re[j];
	R11_im = layers[i].R11_im[j];
	R12_re = layers[i].R12_re[j];
	R12_im = layers[i].R12_im[j];
	R21_re = layers[i].R21_re[j];
	R21_im = layers[i].R21_im[j];
	R22_re = layers[i].R22_re[j];
	R22_im = layers[i].R22_im[j];
	// multiply the transfer matrix of the layer with the
	// R matrix of the layer
// 	M11 = A11*R11 + A12*R21;
// 	M12 = A11*R12 + A12*R22;
// 	M21 = A21*R11 + A11*R21;
// 	M22 = A21*R12 + A11*R22;
	if (complex_A)
	{
	    M11_re = A11_re*R11_re - A11_im*R11_im +
		A12_re*R21_re - A12_im*R21_im;
	    M11_im = A11_re*R11_im + A11_im*R11_re +
		A12_re*R21_im + A12_im*R21_re;
	    M12_re = A11_re*R12_re - A11_im*R12_im +
		A12_re*R22_re - A12_im*R22_im;
	    M12_im = A11_re*R12_im + A11_im*R12_re +
		A12_re*R22_im + A12_im*R22_re;
	    M21_re = A21_re*R11_re - A21_im*R11_im +
		A11_re*R21_re - A11_im*R21_im;
	    M21_im = A21_re*R11_im + A21_im*R11_re +
		A11_re*R21_im + A11_im*R21_re;
	    M22_re = A21_re*R12_re - A21_im*R12_im +
		A11_re*R22_re - A11_im*R22_im;
	    M22_im = A21_re*R12_im + A21_im*R12_re +
		A11_re*R22_im + A11_im*R22_re;
	}
	else
	{
	    M11_re = A11_re*R11_re - A12_im*R21_im;
	    M11_im = A11_re*R11_im + A12_im*R21_re;
	    M12_re = A11_re*R12_re - A12_im*R22_im;
	    M12_im = A11_re*R12_im + A12_im*R22_re;
	    M21_re = - A21_im*R11_im + A11_re*R21_re;
	    M21_im = + A21_im*R11_re + A11_re*R21_im;
	    M22_re = - A21_im*R12_im + A11_re*R22_re;
	    M22_im = + A21_im*R12_re + A11_re*R22_im;
	}
	// store M in R
	R11_re = M11_re; R11_im = M11_im;
	R12_re = M12_re; R12_im = M12_im;
	R21_re = M21_re; R21_im = M21_im;
	R22_re = M22_re; R22_im = M22_im;
	// multiply the L-matrix of the layer with 'R' (which is now AR)
// 	M11 = L11*R11 + L12*R21;
// 	M12 = L11*R12 + L12*R22;
// 	M21 = L21*R11 + L22*R21;
// 	M22 = L21*R12 + L22*R22;
	M11_re = L11_re*R11_re - L11_im*R11_im +
	    + L12_re*R21_re - L12_im*R21_im;
	M11_im = L11_re*R11_im + L11_im*R11_re +
	    L12_re*R21_im + L12_im*R21_re;
	M12_re = L11_re*R12_re - L11_im*R12_im + 
	    L12_re*R22_re - L12_im*R22_im;
	M12_im = L11_re*R12_im + L11_im*R12_re + 
	    L12_re*R22_im + L12_im*R22_re;
	M21_re = L21_re*R11_re - L21_im*R11_im + 
	    L22_re*R21_re - L22_im*R21_im;
	M21_im = L21_re*R11_im + L21_im*R11_re + 
	    L22_re*R21_im + L22_im*R21_re;
	M22_re = L21_re*R12_re - L21_im*R12_im + 
	    L22_re*R22_re - L22_im*R22_im;
	M22_im = L21_re*R12_im + L21_im*R12_re + 
	    L22_re*R22_im + L22_im*R22_re;
	// calculate the optical admittance and the reflectivity
// 	eta = Complex(eta_re[exit_medium_index][j],
// 		      eta_im[exit_medium_index][j]);
	eta_real = eta_re[exit_medium_index][j];
	eta_imag = eta_im[exit_medium_index][j];
//	Y = (M21 + M22*eta) / (M11 + M12*eta);
	Y = Complex(M21_re+M22_re*eta_real-M22_im*eta_imag,
		    M21_im+M22_re*eta_imag+M22_im*eta_real) /
	    Complex(M11_re+M12_re*eta_real-M12_im*eta_imag,
		    M11_im+M12_re*eta_imag+M12_im*eta_real);
	eta0 = Complex(eta_re[incidence_medium_index][j],
		       eta_im[incidence_medium_index][j]);
	reflectivity[j] = (eta0 - Y) / (eta0 + Y);
    }
    // restore the layer's thickness
    layers[i].d = original_thickness;
    UpdateTransferMatrix(i);  // a little waste of time--we could have
			      // stored the original matrix
}

//-------------------------------------------------------------------------

void Coating::Crossover(int crossover_point, const Coating& coating)
{
    int i, j;
    int number_of_layers = NumberOfLayers();

    // some checks
    if (number_of_layers==0) throw "ERROR in Coating::Crossover: \
no layers";
    if (number_of_layers != coating.NumberOfLayers())
    {
	throw "ERROR in Coating::Crossover: the number of layers must \
be the same in both coatings";
    }
    if (NumberOfVariableLayers() != coating.NumberOfVariableLayers())
    {
	throw "ERROR in Coating::Crossover: the number of variable layers \
must be the same in both coatings";
    }
#ifdef DEBUG
    for (i=0; i<number_of_layers; i++)
    {
	if (layers[i].material_name != coating.layers[i].material_name)
	    throw "ERROR in Coating::Crossover: both coatings must be made \
of the same materials";
    }
#endif
    // identify the layer, which actually serves as the crossover point
    for (i=0, j=0; i<number_of_layers && j<crossover_point; i++)
    {
	if (! layers[i].fixed_thickness) j++;
    }
    // it's the i-th layer
    if (i < number_of_layers-1)
    {
	// copy the layers
	for (j=i+1; j<number_of_layers; j++) layers[j] = coating.layers[j];
	// update the L- and R-matrices; it's safer to update the
	// L- matrices from the beginning
	for (j=1; j<number_of_layers; j++) layers[j].UpdateL(layers[j-1]);
	for (j=i; j>=0; j--) layers[j].UpdateR(layers[j+1]);
    }
    Calculate_B_and_C();
    FieldAmplitude(-1, E_0, H_0, eta_0);
}

//-------------------------------------------------------------------------

int Coating::RegisterMaterial(const string& material_name,
			      MaterialRepository& material_repository)
{
    int i, index;
    Real omega, wavelength;
    Real *eta_real, *eta_imag, *pf_real, *pf_imag;
    Complex refractive_index, N0, helper, cos_theta;
    map<string,int>::const_iterator it = material_map.find(material_name);
    if (it==material_map.end())
    {
	index = number_of_materials;
	eta_real = new Real[N];
	eta_imag = new Real[N];
	pf_real = new Real[N];
	pf_imag = new Real[N];
    }
    else
    {
	index = it->second;
	eta_real = eta_re[index];
	eta_imag = eta_im[index];
	pf_real = phase_factor_re[index];
	pf_imag = phase_factor_im[index];
    }

    if (!incidence_medium_name.size())
    {
	throw "ERROR in Coating::RegisterMaterial: \
'incidence_medium_name' is not initialised ";
    }

    // calculate 'eta' and the 'phase_factor'
    for (i=0; i<N; i++)
    {
	omega = omega_min + (omega_max - omega_min) * i/Real(N-1);
	wavelength = 2.0*M_PI*SPEED_OF_LIGHT / omega;
	N0 = material_repository.RefractiveIndex(incidence_medium_name,
						 wavelength);
	if (material_name==incidence_medium_name)
	{
	    refractive_index = N0;
	    cos_theta = cos(angle_of_incidence);
	}
	else
	{
	    refractive_index =
		material_repository.RefractiveIndex(material_name, wavelength);
	    helper = sin(angle_of_incidence) * N0 / refractive_index;
	    cos_theta = cos(ArcSin(helper));
	}
	if (polarisation==TE)
	{
	    helper = refractive_index * cos_theta;
	    eta_real[i] = real(helper);
	    eta_imag[i] = imag(helper);
	}
	else
	{
	    helper = refractive_index / cos_theta;
	    eta_real[i] = real(helper);
	    eta_imag[i] = imag(helper);
	}
	helper = omega * refractive_index * cos_theta / SPEED_OF_LIGHT;
	pf_real[i] = real(helper);
	pf_imag[i] = imag(helper);
    }

    // if it's a new material, insert the arrays into the containers
    if (it==material_map.end())
    {
	eta_re.push_back(eta_real);
	eta_im.push_back(eta_imag);
	phase_factor_re.push_back(pf_real);
	phase_factor_im.push_back(pf_imag);
	material_map.insert(make_pair(material_name, index));
	number_of_materials++;
    }
    return index;
}

//-------------------------------------------------------------------------

void Coating::UpdateTransferMatrix(int i)
{
    int j, k;
    Real delta, helper;
    Complex eta, cdelta, A11, A12, A21, chelper;
    const Complex iu(Real(0), Real(1));
    
    k = material_index[i];
    if (layers[i].complex_A)
    {
	for (j=0; j<N; j++)
	{
	    cdelta = layers[i].d *
		Complex(phase_factor_re[k][j], phase_factor_im[k][j]);
	    A11 = cos(cdelta);
	    chelper = sin(cdelta);
	    eta = Complex(eta_re[k][j], eta_im[k][j]);
	    A12 = iu * chelper / eta;
	    A21 = iu * chelper * eta;
	    layers[i].a11_re[j] = real(A11);
	    layers[i].a11_im[j] = imag(A11);
	    layers[i].a12_re[j] = real(A12);
	    layers[i].a12_im[j] = imag(A12);
	    layers[i].a21_re[j] = real(A21);
	    layers[i].a21_im[j] = imag(A21);
	}
    }
    else
    {
	for (j=0; j<N; j++)
	{
	    delta = layers[i].d * phase_factor_re[k][j];
	    layers[i].a11_re[j] = cos(delta);
	    layers[i].a11_im[j] = 0.0;
	    helper = sin(delta);
	    layers[i].a12_re[j] = 0.0;
	    layers[i].a12_im[j] = helper / eta_re[k][j];
	    layers[i].a21_re[j] = 0.0;
	    layers[i].a21_im[j] = helper * eta_re[k][j];
	}
    }
}

//-------------------------------------------------------------------------

void Coating::RefreshLayers()
{
    int i, j, m, n;

    n = NumberOfLayers();

    // determine for each layer whether the diagonal elements of the
    // transfer matrix are complex (they are complex if the material
    // is absorbing or if we handle the case of the total internal
    // reflection or if the medium of incidence is absorbing)
    for (i=0; i<n; i++)
    {
	layers[i].complex_A = false;
	m = material_index[i];
	for (j=0; j<N; j++)
	{
	    if (abs(phase_factor_im[m][j]) > 1e-7)
	    {
		layers[i].complex_A = true;
		break;
	    }
	}
#ifdef DEBUG
	if (layers[i].complex_A)
	    cout << "complex layer " << i << endl;
#endif
    }

    // calculate the transfer matrices
    for (i=0; i<n; i++) UpdateTransferMatrix(i);

    // calculate the L and R matrices
    for (j=0; j<N; j++)
    {
	layers[0].L11_re[j] = Real(1);
	layers[0].L11_im[j] = 0.0;
	layers[0].L12_re[j] = 0.0;
	layers[0].L12_im[j] = 0.0;
	layers[0].L21_re[j] = 0.0;
	layers[0].L21_im[j] = 0.0;
	layers[0].L22_re[j] = Real(1);
	layers[0].L22_im[j] = 0.0;
    }
    layers[0].complex_L = false;
    layers[0].L_matrix_updated = true;
    for (i=1; i<n; i++) layers[i].UpdateL(layers[i-1]);
    for (j=0; j<N; j++)
    {
	layers[n-1].R11_re[j] = Real(1);
	layers[n-1].R11_im[j] = 0.0;
	layers[n-1].R12_re[j] = 0.0;
	layers[n-1].R12_im[j] = 0.0;
	layers[n-1].R21_re[j] = 0.0;
	layers[n-1].R21_im[j] = 0.0;
	layers[n-1].R22_re[j] = Real(1);
	layers[n-1].R22_im[j] = 0.0;
    }
    layers[n-1].complex_R = false;
    for (i=n-2; i>=0; i--) layers[i].UpdateR(layers[i+1]);
}

//-------------------------------------------------------------------------

void Coating::Calculate_B_and_C()
{
    int j;
    Complex eta; // of the exit medium
    Complex A11, A12, A21;
    Complex R11, R12, R21, R22;
    Complex M11, M12, M21, M22;

    for (j=0; j<N; j++)
    {
	// prepare the complex transfer- and R-matrices of the first layer
	A11 = Complex(layers[0].a11_re[j], layers[0].a11_im[j]);
	A12 = Complex(layers[0].a12_re[j], layers[0].a12_im[j]);
	A21 = Complex(layers[0].a21_re[j], layers[0].a21_im[j]);
	R11 = Complex(layers[0].R11_re[j], layers[0].R11_im[j]);
	R12 = Complex(layers[0].R12_re[j], layers[0].R12_im[j]);
	R21 = Complex(layers[0].R21_re[j], layers[0].R21_im[j]);
	R22 = Complex(layers[0].R22_re[j], layers[0].R22_im[j]);
	// multiply the transfer matrix of the first layer with the
	// R matrix of the layer
	M11 = A11*R11 + A12*R21;
	M12 = A11*R12 + A12*R22;
	M21 = A21*R11 + A11*R21;
	M22 = A21*R12 + A11*R22;
	// calculate B and C
	eta = Complex(eta_re[exit_medium_index][j],
		      eta_im[exit_medium_index][j]);
	B[j] = M11 + M12*eta;
	C[j] = M21 + M22*eta;
    }
}

// void Coating::Calculate_B_and_C()
// {
//     int i, j, n;
//     Complex eta; // of the exit medium
//     Complex A11, A12, A21;
//     Complex m11, m12, m21, m22;
//     Complex M11, M12, M21, M22;
//     n = NumberOfLayers();
//     for (j=0; j<N; j++)
//     {
// 	M11 = Complex(layers[n-1].a11_re[j], layers[n-1].a11_im[j]);
// 	M12 = Complex(layers[n-1].a12_re[j], layers[n-1].a12_im[j]);
// 	M21 = Complex(layers[n-1].a21_re[j], layers[n-1].a21_im[j]);
// 	M22 = M11;
// 	for (i=n-2; i>=0; i--)
// 	{
// 	    A11 = Complex(layers[i].a11_re[j], layers[i].a11_im[j]);
// 	    A12 = Complex(layers[i].a12_re[j], layers[i].a12_im[j]);
// 	    A21 = Complex(layers[i].a21_re[j], layers[i].a21_im[j]);
// 	    m11 = A11*M11 + A12*M21;
// 	    m12 = A11*M12 + A12*M22;
// 	    m21 = A21*M11 + A11*M21;
// 	    m22 = A21*M12 + A11*M22;
// 	    M11 = m11;
// 	    M12 = m12;
// 	    M21 = m21;
// 	    M22 = m22;
// 	}
// 	eta = Complex(eta_re[exit_medium_index][j],
// 		      eta_im[exit_medium_index][j]);
// 	B[j] = M11 + M12*eta;
// 	C[j] = M21 + M22*eta;
//     }
// }

void Coating::FieldAmplitude(int layer, Real z, Complex* E, Complex* H,
			     Complex* _eta)
{
    int j, k;
    int n = NumberOfLayers();
    Complex A11, A12, A21;
    Complex R11, R12, R21, R22;
    Complex M11, M12, M21, M22;
    Complex eta, cdelta, chelper;
    const Complex iu(Real(0), Real(1));

    if (layer<0 || layer>=n)
	throw "ERROR in Coating::FieldAmplitude: layer<0 || layer>=n";
    if (z<Real(0) || z>layers[layer].d)
    {
	throw "\
ERROR in Coating::FieldAmplitude: 'z' beyond the layer boundaries";
    }
    k = material_index[layer];
    layer++; // we calculate the propagation over the distance 'z' separately
    for (j=0; j<N; j++)
    {
	if (layer<n)
	{
	    // prepare the complex transfer- and R-matrices of the layer
	    A11 = Complex(layers[layer].a11_re[j], layers[layer].a11_im[j]);
	    A12 = Complex(layers[layer].a12_re[j], layers[layer].a12_im[j]);
	    A21 = Complex(layers[layer].a21_re[j], layers[layer].a21_im[j]);
	    R11 = Complex(layers[layer].R11_re[j], layers[layer].R11_im[j]);
	    R12 = Complex(layers[layer].R12_re[j], layers[layer].R12_im[j]);
	    R21 = Complex(layers[layer].R21_re[j], layers[layer].R21_im[j]);
	    R22 = Complex(layers[layer].R22_re[j], layers[layer].R22_im[j]);
	    // multiply the transfer matrix of the layer with the
	    // R matrix of the layer
	    M11 = A11*R11 + A12*R21;
	    M12 = A11*R12 + A12*R22;
	    M21 = A21*R11 + A11*R21;
	    M22 = A21*R12 + A11*R22;
	}
	else
	{
	    M11 = 1; M12 = 0; M21 = 0; M22 = 1;
	}
	// take the propagation over the distance 'z' into account
	cdelta = z * Complex(phase_factor_re[k][j], phase_factor_im[k][j]);
	A11 = cos(cdelta);
	chelper = sin(cdelta);
	eta = Complex(eta_re[k][j], eta_im[k][j]);
	if (_eta) _eta[j] = eta;
	A12 = iu * chelper / eta;
	A21 = iu * chelper * eta;
	M11 = A11*M11 + A12*M21;
	M12 = A11*M12 + A12*M22;
	M21 = A21*M11 + A11*M21;
	M22 = A21*M12 + A11*M22;
	// calculate E and H
	eta = Complex(eta_re[exit_medium_index][j],
		      eta_im[exit_medium_index][j]);
	E[j] = M11 + M12*eta;
	H[j] = M21 + M22*eta;
    }
}

//-------------------------------------------------------------------------

void Coating::FieldAmplitude(Real penetration_depth, Complex* E, Complex* H,
			     Complex* _eta)
{
    int i, j, k;
    int n = NumberOfLayers();

    if (penetration_depth<Real(0))
    {
	Complex eta, cdelta, chelper;
	const Complex iu(Real(0), Real(1));
	Complex A11, A12, A21;
	k = incidence_medium_index;
	for (j=0; j<N; j++)
	{
	    cdelta = -penetration_depth *
		Complex(phase_factor_re[k][j], phase_factor_im[k][j]);
	    A11 = cos(cdelta);
	    chelper = sin(cdelta);
	    eta = Complex(eta_re[k][j], eta_im[k][j]);
	    if (_eta) _eta[j] = eta;
	    A12 = iu * chelper / eta;
	    A21 = iu * chelper * eta;
	    E[j] = A11*B[j] + A12*C[j];
	    H[j] = A21*B[j] + A11*C[j];
	}
	return;
    }
    
    if (penetration_depth>StackThickness())
    {
	Complex eta;
	for (j=0; j<N; j++)
	{
	    E[j] = 1;
	    H[j] = Complex(eta_re[exit_medium_index][j],
			   eta_im[exit_medium_index][j]);
	    if (_eta) _eta[j] = H[j];
	}
	return;
    }
    // determine, which layer corresponds to the penetration depth
    // and use the other version of "FieldAmplitude"
    for (i=0; penetration_depth>=Real(0) && i<n; i++)
	penetration_depth -= layers[i].d;
    FieldAmplitude(i-1, -penetration_depth, E, H, _eta);
}

// ---------------------------------------------------------

void Coating::EFieldIntensity(const int N, const int Nx, 
	const Real dx, Real** Z, const bool normalise)
{
    Complex *E, *H, *eta;
    int i, k;
    Real x;

    E = new Complex[N];
    H = new Complex[N];
    eta = new Complex[N];
    
    for (k=0,x=0; k<Nx; k++,x+=dx)
    {
	FieldAmplitude(x, E, H);
	for (i=0; i<N; i++)
	    Z[k][i] = norm(E[i]);
    }

    //normalise the intensity to the intensity of the incident light
    if (normalise)
    {
	for (k=0; k<Nx; k++)
	    for (i=0; i<N; i++)
		Z[k][i] /= norm(E_0[i] + H_0[i]/eta_0[i])/4;
    }

    delete[] E;
    delete[] H;
    delete[] eta;
}


//=========================================================================

void Reflectance(int N, Complex* reflectivity, Real* reflectance)
{
    for (int i = 0; i < N; i++)
    {
	reflectance[i] = norm(reflectivity[i]);
    }
}

//--------------------------------------

void PhaseShift(int N, Complex* reflectivity, Real* phase_shift)
{
    int i;
    Real k; // actually, integer, but I don't care
    Real x;
    // calculate the "raw" phase shift
    for (i=0; i<N; i++)
    {
	phase_shift[i] = arg(reflectivity[i]);
    }

    // correct the first point
    x = phase_shift[1] - phase_shift[0];
    if (x <= 0) k = floor(0.5 - x / M_PI);
    else k = - floor(0.5 + x / M_PI);
    phase_shift[1] += M_PI * k;
    // correct the rest of points
    for (i = 2; i < N; i++)
    {
	x = phase_shift[i] - 2.0*phase_shift[i-1] + phase_shift[i-2];
	if (x <= 0) k = floor(0.5 - x / M_PI);
	else k = - floor(0.5 + x / M_PI);
	phase_shift[i] += M_PI * k;
    }
}

//--------------------------------------

void PhaseToGD(int N, Real omega_step, Real* phase_shift, Real* GD)
{
    Real h2 = 2.0 * omega_step;
    int i;
    for (i = 1; i < N - 1; i++)
    {
	GD[i] = - (phase_shift[i+1] - phase_shift[i-1]) / h2;
    }
    GD[0] = 2.0 * GD[1] - GD[2];
    GD[N-1] = 2.0 * GD[N-2] - GD[N-3];
}

//--------------------------------------

void PhaseToGDD(int N, Real omega_step, Real* phase_shift, Real* GDD)
{
    Real h_sq = square(omega_step);
    int i;
    for (i = 1; i < N - 1; i++)
    {
	GDD[i] =
	    - (phase_shift[i+1] - 2.0*phase_shift[i] + phase_shift[i-1]) / h_sq;
    }
    GDD[0] = 2.0 * GDD[1] - GDD[2];
    GDD[N-1] = 2.0 * GDD[N-2] - GDD[N-3];
}

//=========================================================================

// g++ -DDEBUG -ggdb -Wall coating.cc design.cc material.cc interpolation.cc -o coating
// int main()
// {
//     int i;
//     Complex *reflectivity;
//     Real *R, *GDD, *phase_shift, *wavelengths;
//     try
//     {
// 	Coating coating;
// 	Coating::Parameters parameters;
// 	Design design;
// 	MaterialRepository material_repository("./");
// 	parameters.polarisation = TM;
// 	parameters.angle_of_incidence = 0.0;
// 	parameters.incidence_medium_name = "AIR";
// 	parameters.exit_medium_name = "FUSI";
// 	parameters.N = 500;
// 	parameters.omega_min = 2.0*M_PI*SPEED_OF_LIGHT / 1000.0;
// 	parameters.omega_max = 2.0*M_PI*SPEED_OF_LIGHT / 500.0;
// 	design.ReadFromFile("design1");
// 	coating.SetParameters(parameters, material_repository);
// 	coating.ImportDesign(design, material_repository);
// 	reflectivity = new Complex[parameters.N];
// 	R = new Real[parameters.N];
// 	phase_shift = new Real[parameters.N];
// 	GDD = new Real[parameters.N];
// 	wavelengths = new Real[parameters.N];
// 	coating.Reflectivity(reflectivity);
// 	Reflectance(parameters.N, reflectivity, R);
// 	PhaseShift(parameters.N, reflectivity, phase_shift);
// 	PhaseToGDD(parameters.N, coating.FrequencyStep(),
// 		   phase_shift, GDD);
// 	coating.Wavelengths(wavelengths);
// 	ofstream fs("R.dat");
// 	for (i=0; i<parameters.N; i++)
// 	    fs << wavelengths[i] <<" "<< setprecision(10)
// 	       << setw(20) << R[i]*100.0 << "\n";
// 	fs.close();
// 	fs.open("GDD.dat");
// 	for (i=0; i<parameters.N; i++)
// 	    fs << wavelengths[i] <<" "<< setprecision(10)
// 	       << setw(20) << GDD[i] << "\n";
// 	fs.close();
// 	fs.open("reflectivity.dat");
// 	for (i=0; i<parameters.N; i++)
// 	    fs << wavelengths[i] <<" "<< setprecision(10)
// 	       << setw(20) << real(reflectivity[i])
// 	       <<" "<< imag(reflectivity[i]) << "\n";
// 	fs.close();
// 	coating.ExportDesign(design);
// 	design.SaveToFile("design.xxx");
// 	delete[] reflectivity;
// 	delete[] R;
// 	delete[] phase_shift;
// 	delete[] GDD;
// 	delete[] wavelengths;
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
//     catch(...)
//     {
// 	cerr << "unknown exception" << endl;
// 	return 1;
//     }
//     return 0;
// }
