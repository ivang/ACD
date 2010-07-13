/* AUTHOR: Vladislav S. Yakovlev */

#ifndef __ANALYSIS_HH
#define __ANALYSIS_HH

#include "common.hh"
#include "material.hh"
#include "design.hh"
#include "parameters.hh"

void Analyse(vector<Design>& designs, bool complementary_coatings,
	     const int dots_in_plot,
	     MaterialRepository& material_repository,
	     Parameters& parameters, Variables& variables);

#endif
