// Algorithms for the populational global optimisation

#ifndef __GENETIC_HH
#define __GENETIC_HH

#include "common.hh"

//=========================================================================

class AbstractPopulation
{ /* At first glance it may appear tempting to use something like an
     'AbstractChromosome' in order to provide the global-optimisation
     algorithms with an interface to real problems, and to use one of the
     STL containers to make a population out of the chromosomes However,
     such a solution would have a big problem with the crossover: children
     of a AbstractChromosome would not be allowed to know anything but the
     set of genes of the chromosome they crossover with. Although for some
     applications it may be sufficient, in general you don't want to make
     children with somebody, you only know abstract information about.
     */
public:
    AbstractPopulation() {}
    virtual ~AbstractPopulation() {}
    virtual int Size() = 0; // returns the number of members
    virtual Real MeritFunction(int member) = 0;    
    virtual void KillMember(int member) = 0;
    virtual void Mutate(int member) = 0;
    virtual Real Refine(int member) = 0;
    /* for the post-global optimisation; returns the new value of the
       merit function */

    virtual Real RefineQuickly(int member) = 0;
    /* A fast but incomplete refinement. The new value of the merit
       function is returned. */

    virtual void Sort(void) = 0;
    /* sorts the members in such an order that the best one becomes first and
       the worst one become last */

    virtual int Crossover(int member1, int member2) = 0;
    /* the function creates a new member and adds it to the population;
       the index of the new member is returned */

    virtual bool CheckProximity(int member1, int member2) = 0;
    /* if the two members are detected to be near the same local
       minimum, the function returns true; the function can modify
       (improve) both members, and if they indeed belong to the
       same minimum, they may have the same genes, when the function
       exits */

    virtual void Cut(int n) = 0; // removes 'n' last members of the population
};

//=========================================================================

class PopulationalOptimiser // an abstract class
{
public:
    PopulationalOptimiser();
    virtual ~PopulationalOptimiser();
    virtual void Evolve(AbstractPopulation& population,
			int number_of_generations=1,
			int maximal_population_size=-1) = 0;
    /* "maximal_population_size" is used when the population is allowed to
       grow from generation to generation; if the value of the parameter
       is negative, there is no restriction on the population size */

    void RefineBestMembers(AbstractPopulation& population, int n=1);
    /* at most 'n' best members of the population are placed in
       the beginning of the list and refined; a check is made that all the
       refined members belong to the different local minima, which can reduce
       the size of the population */

    Real ConvergenceFactor(AbstractPopulation& population);

protected:
};

//=========================================================================

class MemeticOptimiser : public PopulationalOptimiser
{
public:
    MemeticOptimiser();
    ~MemeticOptimiser();
    void Evolve(AbstractPopulation& population,
		int number_of_generations=1,
		int maximal_population_size=-1);

protected:
};

//=========================================================================

#endif
