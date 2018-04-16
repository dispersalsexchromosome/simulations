#ifndef INDIVIDU_H_INCLUDED
#define INDIVIDU_H_INCLUDED
#include <bitset>
#include <time.h>
#include "patchs.h"
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <vector>

#include <iostream>
#include <fstream>

struct Individus
{
    // caracs
    bool etat;
    int patch;
    bool dispersant;
    int fecondite;
    int progeniture;
    int sexe_A;
    int sexe_B;
    double disp_A;
    double disp_B;
    std::bitset <32> neutres_A;
    std::bitset <32> neutres_B;
    // fonctions
    void initialisation(const gsl_rng* &gene,int nbr_patch,int chromosexe_a,int chromosexe_b,double disp_init_x,double disp_init_y,int taille_bitset);
    void dispersion(const gsl_rng* &gene,std::vector<Patchs> &parcelles,double &mortalite,int &nbr_patch,int &distance);
    Individus();
    Individus(const gsl_rng* &gene,Individus &papa, Individus &maman,double &taux_de_mut_s,double &importance_des_mut_s,double &taux_de_mut_n,double &sexratio);
};


#endif // INDIVIDU_H_INCLUDED
