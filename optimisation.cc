#include "optimisation.hh"
#include "random.hh"
#include <iostream>

//-------------------------------------------------------------------------

Population::Population(AbstractTarget& target,
		       Real mutation_probability, Real mutation_std) :
    target(target), mutation_probability(mutation_probability),
    mutation_std(mutation_std)
{
}

Population::~Population()
{
}

//-------------------------------------------------------------------------

void Population::Add(Coating& coating)
{
    members.push_back(coating);
    merits.push_back(target.MeritFunction(coating));
}

void Population::Get(int n, Coating& coating)
{
    list<Coating>::iterator it = LocateMember(n);
    coating = *it;
}

//-------------------------------------------------------------------------

void Population::AddRandomDesign(Real minimal_thickness,
				 Real maximal_thickness)
{
    if (Size()==0) throw "ERROR in Population::AddRandomDesign";
    int i, number_of_layers;
    list<Coating>::iterator it, new_member;
    Real* thicknesses;

    // add a copy of the first member to the population
    Coating coating;
    members.push_back(coating);
    new_member = members.end();
    new_member--;
    it = members.begin();
    new_member->operator=(*it);
    // change the variable layers
    number_of_layers = new_member->NumberOfLayers();
    thicknesses = new Real[number_of_layers];
    new_member->GetVariableThicknesses(thicknesses);
    for (i=0; i<number_of_layers; i++)
    {
	thicknesses[i] = UniformRandom(minimal_thickness, maximal_thickness);
    }
    new_member->SetVariableThicknesses(thicknesses);
    merits.push_back(target.MeritFunction(*new_member));
    delete[] thicknesses;
}

//-------------------------------------------------------------------------

Real Population::MeritFunction(int n)
{
#ifdef DEBUG
    if (n<0 || n>=Size()) throw "ERROR in Population::MeritFunction: \
n<0 || n>=Size()";
    list<Coating>::iterator it = LocateMember(n);
    Real merit = target.MeritFunction(*it);
    if (abs(merit-merits[n]) > 1e-5)
    {
	cerr << "A PROBLEM in Population::MeritFunction:\n\
merits[n] = " << merits[n] << "\n\
target.MeritFunction(*it) = " << merit << endl;
	exit(1);
    }
#endif
    return merits[n];
}

//-------------------------------------------------------------------------

void Population::KillMember(int n)
{
    list<Coating>::iterator it_members = LocateMember(n);
    members.erase(it_members);
    vector<Real>::iterator it_merits = merits.begin();
    for (int i=0; i<n && it_merits!=merits.end(); ++it_merits, i++) {}
    if (it_merits==merits.end())
	throw "ERROR in Population::KillMember: 'members' has a size \
different from 'merits'";
    merits.erase(it_merits);
}

//-------------------------------------------------------------------------

void Population::Mutate(int n)
{
    Real x;
    list<Coating>::iterator it = LocateMember(n);
    int number_of_layers = it->NumberOfVariableLayers();
    Real *thicknesses = new Real[number_of_layers];
    it->GetVariableThicknesses(thicknesses);
    for (int i=0; i<number_of_layers; i++)
    {
	x = UniformRandom();
	if (x < mutation_probability)
	{
	    thicknesses[i] = abs(thicknesses[i] +
				 mutation_std * GaussianRandom());
	}
    }
    it->SetVariableThicknesses(thicknesses);
    merits[n] = target.MeritFunction(*it);
    delete[] thicknesses;
}

//-------------------------------------------------------------------------

void Population::Cut(int n)
{ // removes 'n' last members of the population
    for (int i=0; i<n; i++)
    {
	if (Size() > 0)
	{
#ifdef DEBUG
	    cout << "killing the member with the figure of merit "
		 << MeritFunction(Size()-1) << endl;
#endif
	    members.pop_back();
	    merits.pop_back();
	}
    }
}

//-------------------------------------------------------------------------

Real Population::Refine(int n)
{
    list<Coating>::iterator it = LocateMember(n);
    ConjugateGradientOptimisation(*it, target);
    merits[n] = target.MeritFunction(*it);
    return merits[n];
}

//-------------------------------------------------------------------------

Real Population::RefineQuickly(int member)
{   // refine layer by layer using the golden search algorithm
    Real bracketing_step_size = 10.0; // nm
    const int max_steps = 100; // the maximal number of steps, which can
			       // be used to bracket a minimum
    int n, m, counter;
    Real ax, bx, cx;
    Real f1,f2,x0,x1,x2,x3;
    Real prev_merit;
    const Real R = (sqrt(5.0) - 1.0) / 2.0;
    const Real C = 1.0 - R;
    Real merit;
    int* permutation;
    list<Coating>::iterator the_member = LocateMember(member);
    int number_of_layers = the_member->NumberOfVariableLayers();
    Real *thicknesses = new Real[number_of_layers];
    the_member->GetVariableThicknesses(thicknesses);
    merit = target.MeritFunction(*the_member);
    permutation = new int[number_of_layers];
    Permutation(number_of_layers, permutation);
    for (n=0; n<number_of_layers; n++)
    {
	m = permutation[n];
	bx = thicknesses[m];
	// bracket the minimum
	ax = cx = bx;
	prev_merit = merit;
	for (counter=0; counter<max_steps; counter++)
	{
	    ax -= bracketing_step_size;
	    f1 = target.MeritFunction(*the_member, m, ax);
	    if (f1 >= prev_merit) break;
	    prev_merit = f1;
	}
	if (counter==max_steps)
	{   // shouldn't normally happen, but if it happens, just update
	    // the layer thickness
	    thicknesses[m] = ax;
	    the_member->SetVariableThicknesses(thicknesses);
	    merit = prev_merit;
	    continue;
	}
	if (counter > 0)
	{
	    cx = bx;
	    bx = ax + bracketing_step_size;
	    merit = prev_merit;
	}
	else
	{
	    prev_merit = merit;
	    for (counter=0; counter<max_steps; counter++)
	    {
		cx += bracketing_step_size;
		f2 = target.MeritFunction(*the_member, m, cx);
		if (f2 >= prev_merit) break;
		prev_merit = f2;
	    }
	    if (counter==max_steps)
	    {   // shouldn't normally happen, but if it happens, just update
		// the layer thickness
		thicknesses[m] = cx;
		the_member->SetVariableThicknesses(thicknesses);
		merit = prev_merit;
		continue;
	    }
	    if (f1 == merit && f2 == merit)
	    {	// the merit function seems to be independent on that variable
		continue;
	    }
	}
	// perform the golden section search
	x0=ax;
	x3=cx;
	if (abs(cx-bx) > abs(bx-ax))
	{
	    x1=bx;
	    x2=bx+C*(cx-bx);
	    f1 = merit;
	    f2 = target.MeritFunction(*the_member, m, x2);
	}
	else
	{
	    x2=bx;
	    x1=bx-C*(bx-ax);
	    f2 = merit;
	    f1 = target.MeritFunction(*the_member, m, x1);
	}
	while (abs(x3-x0) > Real(1))
	{
	    if (f2 < f1)
	    {
		x0 = x1;
		x1 = x2;
		x2 = R*x1+C*x3;
		f1 = f2;
		f2 = target.MeritFunction(*the_member, m, x2);
	    }
	    else
	    {
		x3 = x2;
		x2 = x1;
		x1 = R*x2+C*x0;
		f2 = f1;
		f1 = target.MeritFunction(*the_member, m, x1);
	    }
	}
	if (f1 < f2)
	{
	    thicknesses[m] = x1;
	    the_member->SetVariableThicknesses(thicknesses);	    
	    merit = f1;
	}
	else
	{
	    thicknesses[m] = x2;
	    the_member->SetVariableThicknesses(thicknesses);
	    merit = f2;
	}
    } // for (n = 0; n < number_of_layers; n++)
    merits[member] = target.MeritFunction(*the_member);
    delete[] permutation;
    delete[] thicknesses;
    return merits[member];
}

//-------------------------------------------------------------------------

void Population::Sort(void)
{
    int i, N;
    Real x;
    bool swapped;
    list<Coating>::iterator it, it1, it2;
    N = Size();
    if (N <= 1) return;
    // double-pass bubble sort
    do
    {
	swapped = false;
	it1=members.begin();
	it2 = it1;
	++it2;
	for (i=0; i<N-1; i++)
	{
	    if (merits[i+1] < merits[i])
	    {
		x = merits[i]; merits[i] = merits[i+1]; merits[i+1] = x;
		members.splice(it1, members, it2);
		it = it1; it1 = it2; it2 = it;
		swapped = true;
	    }
	    it1 = it2++;
	}
	if (! swapped) break;
	it1=members.end();
	--it1;
	it2 = it1;
	--it2;
	for (i=N-1; i>0; i--)
	{
	    if (merits[i-1] > merits[i])
	    {
		x = merits[i]; merits[i] = merits[i-1]; merits[i-1] = x;
		members.splice(it2, members, it1);
		it = it1; it1 = it2; it2 = it;
		swapped = true;
	    }
	    it1 = it2--;
	}
    } while(swapped);
#ifdef DEBUG
    swapped = false;
    for (i=0; i<N-1; i++)
    {
	if (merits[i] > merits[i+1])
	{
	    swapped = true;
	    break;
	}
    }
    if (swapped)
    {
	for (i=0; i<N; i++)
	{
	    cout << i <<" "<< merits[i] << "\n";
	}
	exit(1);
    }
#endif
}

//-------------------------------------------------------------------------

int Population::Crossover(int member1, int member2)
{
    int crossover_point, number_of_layers, n;
    list<Coating>::iterator it, new_member;
    // add a copy of the first member to the population
    Coating coating;
    members.push_back(coating);
    merits.push_back(Real(0));
    new_member = members.end();
    new_member--;
    it = LocateMember(member1);
    new_member->operator=(*it);
    // determine the crossover point
    number_of_layers = new_member->NumberOfLayers();
    if (number_of_layers > 2)
    {
	crossover_point = RandomInteger(1, number_of_layers-2);
#ifdef DEBUG
	cout << "[" << crossover_point << "] ";
#endif
    }
    else crossover_point = number_of_layers-1;
    // crossover;
    it = LocateMember(member2);
    new_member->Crossover(crossover_point, *it);
    n = Size()-1;
    merits[n] = target.MeritFunction(*new_member);
    return n;
}

//-------------------------------------------------------------------------

bool Population::CheckProximity(int member1, int member2)
{
    const Real threshold = 1e-4; // if the reflectivities differ by less
				 // than this value, the designs are
				 // considered to be identical
    int j, N;
    Real discrepancy;
    bool designs_are_identical;
    list<Coating>::iterator iterator1 = LocateMember(member1);
    list<Coating>::iterator iterator2 = LocateMember(member2);
    Complex *reflectivity1, *reflectivity2;

    N = iterator1->NumberOfFrequencies();
    if (N != iterator2->NumberOfFrequencies())
    {
	throw "ERROR in Population::CheckProximity: \n\
the designs have different numbers of analysis frequencies";
    }
    reflectivity1 = new Complex[N];
    reflectivity2 = new Complex[N];
    // compare the reflectivities of the coatings
    iterator1->Reflectivity(reflectivity1);
    iterator2->Reflectivity(reflectivity2);
    designs_are_identical = true;
    for (j=0; j<N; j++)
    {
	discrepancy = norm(reflectivity2[j] - reflectivity1[j]);
	if (discrepancy > threshold)
	{
	    designs_are_identical = false;
	    break;
	}
    }
    if (designs_are_identical)
    {
	// copy the design with a smaller figure of merit to the design
	// with a larger figure of merit
	if (merits[member1] < merits[member2])
	{
	    *iterator2 = *iterator1;
	    merits[member2] = merits[member1];
	}
	else
	{
	    *iterator1 = *iterator2;
	    merits[member1] = merits[member2];
	}
    }
	
    delete[] reflectivity1;
    delete[] reflectivity2;
    return designs_are_identical;
    /* quite some time is wasted in this function, because the reflectivity
       of the same design is calculated many times; however, the time wasted
       is very small in comparison with the time spent in the optimisation
       subroutines
    */
}

//-------------------------------------------------------------------------

list<Coating>::iterator Population::LocateMember(int n)
{
    int i;
    list<Coating>::iterator it = members.begin();
    for (i=0; i<n && it!=members.end(); ++it, i++) {}
    if (it==members.end())
    {
	throw "ERROR in Population::LocateMember";
    }
    return it;
}

//-------------------------------------------------------------------------

Real Population::BestMerit(int* member)
{
    int i, m, n;
    Real best_merit;
    n = Size();
    if (n<=0) throw "ERROR in Population::BestMerit";
    m = 0;
    best_merit = merits[0];
    for (i=1; i<n; i++)
    {
	if (best_merit > merits[i])
	{
	    best_merit = merits[i];
	    m = i;
	}
    }
    if (member) *member = m;
    return best_merit;
}

//=========================================================================
