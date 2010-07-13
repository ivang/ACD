// all you need to optimise coatings

#ifndef __OPTIMISATION_HH
#define __OPTIMISATION_HH

#include "common.hh"
#include "coating.hh"
#include "genetic.hh"
#include "local.hh"
#include "target.hh"
#include <list>


//-------------------------------------------------------------------------
// The implementation of the population abstraction for the
// global optimisation of coatings

class Population: public AbstractPopulation
{
public:
    Population(AbstractTarget& target, Real mutation_probability,
	       Real mutation_std);
    /* 'mutation_probability' is the probability that a particular
       gene will be mutated, when "void Mutate(int)" is called;
       in this implementation the value of the thickness perturbation
       is normally distributed with the standard deviation "mutation_std" */

    ~Population();

    void Add(Coating& coating);
    // makes a copy of the coating and stores it in the populations

    void AddRandomDesign(Real minimal_thickness=50.0,
			 Real maximal_thickness=300.0);

    void Get(int n, Coating& coating);
    // copies the n-th coating to the given object

    int Size() { return members.size(); }
    Real MeritFunction(int member);
    void KillMember(int member);
    void Mutate(int member);
    Real Refine(int member); // for the post-global optimisation

    Real RefineQuickly(int member);
    /* A fast but incomplete refinement. The new value of the merit
       function is returned. */

    void Sort(void);
    /* sorts the members in such an order that the best one becomes first and
       the worst one become last */

    int Crossover(int member1, int member2);
    /* the function creates a new member and adds it to the population;
       the index of the new member is returned */

    bool CheckProximity(int member1, int member2);
    /* if the two members are detected to be near the same local
       minimum, the function returns true; the function can modify
       (improve) both members, and if they indeed belong to the
       same minimum, they may have the same genes, when the function
       exits */

    void Cut(int n); // removes 'n' last members of the population

    Real BestMerit(int* member=NULL);

private:
    AbstractTarget& target;
    list<Coating> members;
    vector<Real> merits;
    list<Coating>::iterator LocateMember(int n);
    Real mutation_probability;
    Real mutation_std;
};

#endif
