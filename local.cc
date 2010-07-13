#include "local.hh"
#include "gsl/gsl_multimin.h"

extern "C"
{
#	include "numrec.h"
}

//-------------------------------------------------------------------------
// the following variables are initialised when a local optimisation
// starts

AbstractTarget* global_target = NULL;
Coating* optimised_coating = NULL;
vector<Coating>* optimised_set_of_coatings = NULL;

//-------------------------------------------------------------------------
// the following functions are necessary for the external
// optimisation functions

void UpdateOptimisedCoatings(double* optimisation_parameters)
{ // calls 'SetVariableThicknesses' for each of the coatings
    int i, j, k, n, N;
    Coating* c;
    Real* thicknesses;
    n = optimised_set_of_coatings->size();
    j = 0;
    for (i=0; i<n; i++)
    {
	c = &((*optimised_set_of_coatings)[i]);
	N = c->NumberOfVariableLayers();
	thicknesses = new Real[N];
	for (k=0; k<N; k++)
	    thicknesses[k] = abs(optimisation_parameters[j++]);
	c->SetVariableThicknesses(thicknesses);
	delete[] thicknesses;
    }
}


double funk1(const gsl_vector *arguments, void *params)
{

    if (optimised_coating)
    { // we assume that a single mirror is being optimised
	int n = optimised_coating->NumberOfVariableLayers();
	Real* thicknesses = new Real[n];
	for (int i=0; i<n; i++) thicknesses[i] = abs(gsl_vector_get(arguments, i));
	optimised_coating->SetVariableThicknesses(thicknesses);
	delete[] thicknesses;
	return global_target->MeritFunction(*optimised_coating);
    }
    else if (optimised_set_of_coatings)
    { // simultaneous optimisation of a set of complementary mirrors
	int i, number_of_coatings = optimised_set_of_coatings->size();
    	double args[number_of_coatings];

	for (i=0; i<number_of_coatings; i++) args[i] = gsl_vector_get(arguments, i);

	UpdateOptimisedCoatings(&args[1]);
        return global_target->MeritFunction(*optimised_set_of_coatings);
    }
    else throw "ERROR in 'funk1'";
    return 0.0;
}


void dfunk(const gsl_vector *arguments, void *params, gsl_vector *gradient)
{
    int i;
    int n;
    Real* real_array;
    if (optimised_coating)
    { // we assume that a single mirror is being optimised
	n = optimised_coating->NumberOfVariableLayers();
	real_array = new Real[n];
	for (i=0; i<n; i++) real_array[i] = abs(gsl_vector_get(arguments, i));
	optimised_coating->SetVariableThicknesses(real_array);
	global_target->MeritGradient(*optimised_coating, real_array);
    }
    else if (optimised_set_of_coatings)
    { // simultaneous optimisation of a set of complementary mirrors
	int number_of_coatings = optimised_set_of_coatings->size();
	double args[number_of_coatings];

	n = 0;
	for (i=0; i<number_of_coatings; i++) {
	    n += (*optimised_set_of_coatings)[i].NumberOfVariableLayers();
	    args[i] = gsl_vector_get(arguments, i);
	}
	real_array = new Real[n];
	UpdateOptimisedCoatings(&args[1]);
        global_target->MeritGradient(*optimised_set_of_coatings, real_array);
    }
    else throw "ERROR in 'dfunk'";
    for (i=0; i<n; i++) gsl_vector_set(gradient, i, real_array[i]);
    delete[] real_array;
}

void fdfunk(const gsl_vector *arguments, void *params, double *f, gsl_vector *gradient)
{
	*f = funk1(arguments, params);
	dfunk(arguments, params, gradient);
}

//-------------------------------------------------------------------------

void ConjugateGradientOptimisation(Coating& coating,
				   AbstractTarget& target)
{
    const int params = NULL;    // GSL requires such parameters. Not needed
   				// in our program, though.

    const gsl_multimin_fdfminimizer_type *T; 
    gsl_multimin_fdfminimizer *gradient_minimizer = NULL;
    gsl_vector *x;
    gsl_multimin_function_fdf merit_min;

    int status;
    int i, N;
    const float ftol = 1e-8; // the fractional tolerance
    const float step = 0.01; // initial step size
    
    size_t iter = 0; // the number of iterations that were performed
    
    N = coating.NumberOfVariableLayers();
    if (N<=0) return;
    
    global_target = &target;
    optimised_coating = &coating;

    // initialise the starting point
    Real* thicknesses = new Real[N];
    coating.GetVariableThicknesses(thicknesses);

    merit_min.n = N;
    merit_min.f = funk1;
    merit_min.df = dfunk;
    merit_min.fdf = fdfunk;
    merit_min.params = params;

    x = gsl_vector_alloc(N);
    for (i=0; i<N; i++) gsl_vector_set(x, i, thicknesses[i]);
    
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    gradient_minimizer = gsl_multimin_fdfminimizer_alloc(T, N);

    // Optimize
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

    // Set optimized thicknesses
    for (i=0; i<N; i++) 
	thicknesses[i] = abs(gsl_vector_get(gradient_minimizer->x, i));
    coating.SetVariableThicknesses(thicknesses);
    
    delete[] thicknesses;
    global_target = NULL;
    optimised_coating = NULL;
}

void ConjugateGradientOptimisation(vector<Coating>& coatings,
				   AbstractTarget& target)
{
    const int params = NULL;    // GSL requires such parameters. Not needed
   				// in our program, though.

    const gsl_multimin_fdfminimizer_type *T; 
    gsl_multimin_fdfminimizer *gradient_minimizer = NULL;
    gsl_vector *x;
    gsl_multimin_function_fdf merit_min;

    int status;
    int i, j, k, n, m, N;
    const float ftol = 1e-8; // the fractional tolerance
    const float step = 0.01; // initial step size
    
    size_t iter = 0; // the number of iterations that were performed
    Real* thicknesses;
    n = coatings.size();
    for (i=0, N=0; i<n; i++) N += coatings[i].NumberOfVariableLayers();
    
    if (N<=0) return;
    
    global_target = &target;
    optimised_set_of_coatings = &coatings;

    // initialise the starting point
    merit_min.n = N;
    merit_min.f = funk1;
    merit_min.df = dfunk;
    merit_min.fdf = fdfunk;
    merit_min.params = params;

    x = gsl_vector_alloc(N);
    j = 0;
    for (i=0; i<n; i++)
    {
	m = coatings[i].NumberOfVariableLayers();
	thicknesses = new Real[m];
	coatings[i].GetVariableThicknesses(thicknesses);
	for (k=0; k<m; k++) gsl_vector_set(x, j++, thicknesses[k]);
	delete[] thicknesses;
    }
    
    T = gsl_multimin_fdfminimizer_conjugate_fr;
    gradient_minimizer = gsl_multimin_fdfminimizer_alloc(T, N);

    // Optimize
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

    // Set optimized thicknesses
    for (i=0; i<N; i++) 
	thicknesses[i] = abs(gsl_vector_get(gradient_minimizer->x, i));

    UpdateOptimisedCoatings(&thicknesses[1]);
    
    delete[] thicknesses;
    global_target = NULL;
    optimised_coating = NULL;
}

//-------------------------------------------------------------------------

void SimplexOptimisation(Coating& coating, AbstractTarget& target, Real step)
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
    int i, N;

    const float ftol = 1e-4; // the fractional tolerance
    N = coating.NumberOfVariableLayers();
    if (N<=0) return;

    global_target = &target;
    optimised_coating = &coating;
    
    // Get the initial thicknesses 
    Real* thicknesses = new Real[N];
    coating.GetVariableThicknesses(thicknesses);

    // Starting point
    x = gsl_vector_alloc(N);
    for (i=0; i<N; i++) gsl_vector_set(x, i, thicknesses[i]); 

    steps = gsl_vector_alloc(N);
    gsl_vector_set_all(steps, step);
	
    // Initialize
    merit_min.n = N;
    merit_min.f = funk1;
    merit_min.params = params;

    simplex = gsl_multimin_fminimizer_alloc(T, N);
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

    } while (status == GSL_CONTINUE && iter < 10000);

    for (i=0; i<N; i++) 
	thicknesses[i] = abs(gsl_vector_get(simplex->x, i));
    coating.SetVariableThicknesses(thicknesses);

    delete[] thicknesses;
    global_target = NULL;
    optimised_coating = NULL;
}


void SimplexOptimisation(vector<Coating>& coatings,
	AbstractTarget& target, Real step)
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
    int i, j, k, m, n, N;

    const float ftol = 1e-4; // the fractional tolerance
    n = coatings.size();
    for (i=0, N=0; i<n; i++) N += coatings[i].NumberOfVariableLayers();
    
    if (N<=0) return;

    global_target = &target;
    optimised_set_of_coatings = &coatings;
    Real* thicknesses = new Real[N];

    // Starting point
    x = gsl_vector_alloc(N);
    for (i=0; i<N; i++) gsl_vector_set(x, i, thicknesses[i]); 

    steps = gsl_vector_alloc(N);
    gsl_vector_set_all(steps, step);
	
    // Initialize
    merit_min.n = N;
    merit_min.f = funk1;
    merit_min.params = params;

    x = gsl_vector_alloc(N);
    j = 0;
    for (i=0; i<n; i++)
    {
	m = coatings[i].NumberOfVariableLayers();
	thicknesses = new Real[m];
	coatings[i].GetVariableThicknesses(thicknesses);
	for (k=0; k<m; k++) gsl_vector_set(x, j++, thicknesses[k]);
	delete[] thicknesses;
    }
    
    simplex = gsl_multimin_fminimizer_alloc(T, N);
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

    } while (status == GSL_CONTINUE && iter < 10000);

    for (i=0; i<N; i++) 
	thicknesses[i] = abs(gsl_vector_get(simplex->x, i));

    UpdateOptimisedCoatings(&thicknesses[1]);
    
    delete[] thicknesses;
    global_target = NULL;
    optimised_coating = NULL;
}
