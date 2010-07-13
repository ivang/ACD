#include "genetic.hh"
#include "random.hh"
#include <cstdlib>

//=========================================================================

PopulationalOptimiser::PopulationalOptimiser()
{
}

PopulationalOptimiser::~PopulationalOptimiser() {}

//-------------------------------------------------------------------------

void PopulationalOptimiser::
RefineBestMembers(AbstractPopulation& population, int n)
{
    int i, j, N;
    bool member_killed;
    // refine the best member
    population.Sort();
    population.Refine(0);
    // for each of the following members make sure that they don't
    // belong to the same local minimum with one of the previously
    // optimised members
    N = population.Size();
    for (i=1; i<n && i<N; i++)
    {
	do
	{
	    population.Refine(i);
	    member_killed = false;
	    for (j=0; j<i; j++)
	    {
		if (population.CheckProximity(i, j))
		{
		    population.KillMember(i);
		    member_killed = true;
		    N--;
		    break;
		}
	    }
	} while (member_killed && i<N);
    }
}

//-------------------------------------------------------------------------

Real PopulationalOptimiser::
ConvergenceFactor(AbstractPopulation& population)
{
    Real f, f_min, f_max;
    int i;
    int N = population.Size();
    if (N < 1) return 0.0;
    // determine the largest and the smallest value of the merit
    // function in the population
    f_min = f_max = population.MeritFunction(0);
    for (i=1; i<N; i++)
    {
	f = population.MeritFunction(i);
	if (f_min > f) f_min = f;
	if (f_max < f) f_max = f;
    }
    // return the convergence factor
    if (f_max==Real(0) && f_min==Real(0)) return 0.0;
    return (f_max - f_min) / (f_max + f_min);
}

//=========================================================================

MemeticOptimiser::MemeticOptimiser() : PopulationalOptimiser()
{
}

MemeticOptimiser::~MemeticOptimiser()
{
}

//-------------------------------------------------------------------------

void MemeticOptimiser::Evolve(AbstractPopulation& population,
			      int number_of_generations,
			      int maximal_population_size)
{
    int i, generation, sum_of_ranks, x;
    int member1, member2;
    int population_size;
#ifdef DEBUG
    Real merit;
#endif

    for (generation=0; generation<number_of_generations; generation++)
    {
	population.Sort();
	population_size = population.Size();
	// reduce the size of the population if necessary
	if (population_size > maximal_population_size)
	{
	    population.Cut(population_size - maximal_population_size);
	    population_size = population.Size();
	}
	if (population_size < 2) return;
	// choose two members for the crossover using the rank selection
	sum_of_ranks = population_size * (population_size + 1) / 2;
	x = RandomInteger(0, sum_of_ranks-1);
	for (i=0; i<population_size; i++)
	{
	    x -= (population_size - i);
	    if (x <= 0) break;
	}
	member1 = i;
	do
	{
	    x = RandomInteger(0, sum_of_ranks-1);
	    for (i=0; i<population_size; i++)
	    {
		x -= (population_size - i);
		if (x <= 0) break;
	    }
	    member2 = i;
	} while (member1 == member2);
#ifdef DEBUG
	merit = population.MeritFunction(member1);
	cout <<	merit << ",";
	merit = population.MeritFunction(member2);
	cout <<	merit << " (" << member1 << "," << member2
	     << ") --> " << flush;
#endif
	// crossover
	i = population.Crossover(member1, member2);
#ifdef DEBUG
	merit = population.MeritFunction(i);
	cout <<	merit << " --> ";
#endif
	// mutation
	population.Mutate(i);
#ifdef DEBUG
	merit = population.MeritFunction(i);
	cout <<	merit << " --> ";
#endif
	// partial refinement
#ifdef DEBUG
	merit = population.RefineQuickly(i);
	cout <<	merit << endl;
	cout << "checking the merit function: ";
	cout << abs(merit - population.MeritFunction(i)) << endl;
#else
	population.RefineQuickly(i);
#endif
    }
}

//=========================================================================


/* TESTING BLOCK

// g++ -Wall -ggdb -DDEBUG -o genetic genetic.cc

#include <vector>
#include <list>

class Population: public AbstractPopulation
{
private:
    int number_of_genes;
    list< vector<Real> > members;
    list< vector<Real> >::iterator LocateMember(int n);
public:
    Population(int number_of_genes) : number_of_genes(number_of_genes) {}
    ~Population() {}
    int Size() { return members.size(); }
    Real Refine(int member) { return MeritFunction(member); }
    Real RefineQuickly(int member) { return MeritFunction(member); }
    bool CheckProximity(int member1, int member2) {return false;}
    int AddMember(const vector<Real>& gene);
    int GenerateRandomMember();
    Real MeritFunction(int member);
    Real MeritFunction(const list< vector<Real> >::iterator&);
    void KillMember(int member);
    void Mutate(int member);
    void Sort(void);
    int Crossover(int member1, int member2);
    void Cut(int n);
};


int Population::AddMember(const vector<Real>& genes)
{
    if ((int)genes.size() != number_of_genes)
	throw "gene.size() != number_of_genes";
    members.push_back(genes);
    return Size()-1;
}

list< vector<Real> >::iterator Population::LocateMember(int n)
{
    int i;
    list< vector<Real> >::iterator it = members.begin();
    for (i=0; i<n && it!=members.end(); ++it, i++) {}
    return it;
}

Real Population::MeritFunction(int n)
{
    list< vector<Real> >::iterator it = LocateMember(n);
    return MeritFunction(it);
}

Real Population::MeritFunction(const list< vector<Real> >::iterator& it)
{
    Real result = 0.0;
    if (it == members.end()) return 0.0;
    for (int i=0; i<number_of_genes; i++) result += square((*it)[i]);
    return result;
}

void Population::KillMember(int n)
{
    list< vector<Real> >::iterator it = LocateMember(n);
    members.erase(it);
}

void Population::Mutate(int n)
{
    Real x;
    Real probability = 0.01;
    list< vector<Real> >::iterator it = LocateMember(n);
    if (it == members.end()) return;
    for (int i=0; i<number_of_genes; i++)
    {
	x = UniformRandom(0, 1)
	if (x > probability)
	    (*it)[i] = UniformRandom(-10, 10);
    }
}

void Population::Sort(void)
{
    int i, N;
    Real x;
    bool swapped;
    list< vector<Real> >::iterator it, it1, it2;
    // make a list of merit functions
    N = Size();
    if (N <= 1) return;
    vector<Real> merits(N);
    for (i=0, it=members.begin(); i<N; ++it, i++)
	merits[i] = MeritFunction(it);
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
}

int Population::Crossover(int member1, int member2)
{
    int i, crossover_point;
    vector<Real> new_genes(number_of_genes);
    list< vector<Real> >::iterator it1, it2;
    it1 = LocateMember(member1);
    it2 = LocateMember(member2);
    // choose the crossover point
    crossover_point = RandomInteger(0, number_of_genes-1);
    // copy the genes
    for (i=0; i<=crossover_point; i++) new_genes[i] = (*it1)[i];
    for (i=crossover_point+1; i<number_of_genes; i++)
	new_genes[i] = (*it1)[i];
    // add the new member
    return AddMember(new_genes);
}

void Population::Cut(int n)
{
    list< vector<Real> >::iterator it = members.end();
    for (int i=0; i<n; i++) --it;
    members.erase(it, members.end());
}

int Population::GenerateRandomMember()
{
    vector<Real> genes(number_of_genes);
    for (int i=0; i<number_of_genes; i++)
    {
	genes[i] = UniformRandom(-10, 10);
    }
    return AddMember(genes);
}

int main()
{
    const int number_of_genes = 10;
    const int population_size = 10;
    const int number_of_generations = 1000;
    int i, j;
    try
    {
	Population population(number_of_genes);
	// create the initial population
	cout << "initial merits:" << endl;
	for (i=0; i<population_size; i++)
	{
	    j = population.GenerateRandomMember();
	    cout << population.MeritFunction(j) << endl;
	}
	// optimise
	MemeticOptimiser optimiser;
	optimiser.Evolve(population, number_of_generations, population_size);
	cout << "merits after " << number_of_generations
	     << " generations" << endl;
	for (i=0; i<population_size; i++)
	{
	    cout << population.MeritFunction(i) << endl;
	}
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

*/
