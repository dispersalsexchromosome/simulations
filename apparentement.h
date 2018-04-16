#ifndef APPARENTEMENT_H_INCLUDED
#define APPARENTEMENT_H_INCLUDED
#include "Individu.h"
#include <bitset>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

double apparentement(const gsl_rng* &gene, Individus &focal, Individus &compare);

#endif // APPARENTEMENT_H_INCLUDED
