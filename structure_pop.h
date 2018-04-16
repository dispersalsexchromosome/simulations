#ifndef STRUCTURE_POP_H_INCLUDED
#define STRUCTURE_POP_H_INCLUDED
#include "Individu.h"
#include "apparentement.h"
#include <vector>
#include <iostream>
#include <fstream>

void structure_pop_avd(const gsl_rng* &gene,std::vector<Individus> &males,std::vector<Individus> &femelles,char *nomsortie);
void structure_pop_apd(const gsl_rng* &gene,std::vector<Individus> &males,std::vector<Individus> &femelles,char *nomsortie);

#endif // STRUCTURE_POP_H_INCLUDED
