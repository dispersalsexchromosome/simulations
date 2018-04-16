#include "Individu.h"
#include "apparentement.h"
#include <vector>
#include <iostream>
#include <fstream>

void structure_pop_avd(const gsl_rng* &gene,std::vector<Individus> &males,std::vector<Individus> &femelles,char *nomsortie) {

    std::ofstream sortie_apparentement(nomsortie,std::ios::app);

    const int nbr_comparaison = 1000;
    int nbr_males = males.size();
    int nbr_femelles = femelles.size();
    int a=0;
    int b=0;
    double appa_mm_g = 0;
    double appa_mf_g = 0;
    double appa_ff_g = 0;
    double appa_mm_l = 0;
    double appa_mf_l = 0;
    double appa_ff_l = 0;
    for (int n=0;n<nbr_comparaison;n++) {
        //INTER PATCH
        // apparentement male male inter patch
        a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_males);
        //while (males[a].patch==males[b].patch) {
        //    b = gsl_rng_uniform_int(gene,nbr_males);
        // }
        appa_mm_g = appa_mm_g+apparentement(gene,males[a],males[b]);
        // apparentement male femelle (et vice verca) inter patch
        a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_femelles);
        //while (males[a].patch==femelles[b].patch) {
        //    b = gsl_rng_uniform_int(gene,nbr_femelles);
        //}
        appa_mf_g = appa_mf_g+apparentement(gene,males[a],femelles[b]);
        // apparentement femelle femelle inter patch
        a = gsl_rng_uniform_int(gene,nbr_femelles);b = gsl_rng_uniform_int(gene,nbr_femelles);
        //while (femelles[a].patch==femelles[b].patch) {
        //    b = gsl_rng_uniform_int(gene,nbr_femelles);
        //}
        appa_ff_g = appa_ff_g+apparentement(gene,femelles[a],femelles[b]);


        //INTRA PATCH
        // apparentement male male intra patch
        a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_males);
        while (males[a].patch!=males[b].patch || a==b) {
            a = gsl_rng_uniform_int(gene,nbr_males) ; b = gsl_rng_uniform_int(gene,nbr_males);
        }
        appa_mm_l = appa_mm_l+apparentement(gene,males[a],males[b]);

        // apparentement male femelle (et vice verca) intra patch
        a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_femelles);
        while (males[a].patch!=femelles[b].patch || a==b) {
            a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_femelles);
        }
        appa_mf_l = appa_mf_l+apparentement(gene,males[a],femelles[b]);
        // apparentement femelle femelle intra patch
        a = gsl_rng_uniform_int(gene,nbr_femelles);b = gsl_rng_uniform_int(gene,nbr_femelles);
        while (femelles[a].patch!=femelles[b].patch || a==b) {
            a = gsl_rng_uniform_int(gene,nbr_femelles);b = gsl_rng_uniform_int(gene,nbr_femelles);
        }
        appa_ff_l = appa_ff_l+apparentement(gene,femelles[a],femelles[b]);

    }
    appa_mm_g = appa_mm_g/nbr_comparaison;
    appa_mf_g = appa_mf_g/nbr_comparaison;
    appa_ff_g = appa_ff_g/nbr_comparaison;
    appa_mm_l = appa_mm_l/nbr_comparaison;
    appa_mf_l = appa_mf_l/nbr_comparaison;
    appa_ff_l = appa_ff_l/nbr_comparaison;

    sortie_apparentement << appa_mm_g <<"\t"<< appa_mf_g <<"\t"<< appa_ff_g <<"\t"<< appa_mm_l <<"\t"<< appa_mf_l <<"\t"<< appa_ff_l <<"\t";
    sortie_apparentement.close();
}


void structure_pop_apd(const gsl_rng* &gene,std::vector<Individus> &males,std::vector<Individus> &femelles,char *nomsortie) {

    std::ofstream sortie_apparentement(nomsortie,std::ios::app);

    const int nbr_comparaison = 1000;
    int nbr_males = males.size();
    int nbr_femelles = femelles.size();
    int a=0;
    int b=0;
    double appa_mm_g = 0;
    double appa_mf_g = 0;
    double appa_ff_g = 0;
    double appa_mm_l = 0;
    double appa_mf_l = 0;
    double appa_ff_l = 0;
    for (int n=0;n<nbr_comparaison;n++) {
        //INTER PATCH
        // apparentement male male inter patch
        a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_males);
        //while (males[a].patch==males[b].patch) {
        //    b = gsl_rng_uniform_int(gene,nbr_males);
        // }
        appa_mm_g = appa_mm_g+apparentement(gene,males[a],males[b]);
        // apparentement male femelle (et vice verca) inter patch
        a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_femelles);
        //while (males[a].patch==femelles[b].patch) {
        //    b = gsl_rng_uniform_int(gene,nbr_femelles);
        //}
        appa_mf_g = appa_mf_g+apparentement(gene,males[a],femelles[b]);
        // apparentement femelle femelle inter patch
        a = gsl_rng_uniform_int(gene,nbr_femelles);b = gsl_rng_uniform_int(gene,nbr_femelles);
        //while (femelles[a].patch==femelles[b].patch) {
        //    b = gsl_rng_uniform_int(gene,nbr_femelles);
        //}
        appa_ff_g = appa_ff_g+apparentement(gene,femelles[a],femelles[b]);


        //INTRA PATCH
        // apparentement male male intra patch
        a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_males);
        while (males[a].patch!=males[b].patch || a==b) {
            a = gsl_rng_uniform_int(gene,nbr_males) ; b = gsl_rng_uniform_int(gene,nbr_males);
        }
        appa_mm_l = appa_mm_l+apparentement(gene,males[a],males[b]);

        // apparentement male femelle (et vice verca) intra patch
        a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_femelles);
        while (males[a].patch!=femelles[b].patch || a==b) {
            a = gsl_rng_uniform_int(gene,nbr_males);b = gsl_rng_uniform_int(gene,nbr_femelles);
        }
        appa_mf_l = appa_mf_l+apparentement(gene,males[a],femelles[b]);
        // apparentement femelle femelle intra patch
        a = gsl_rng_uniform_int(gene,nbr_femelles);b = gsl_rng_uniform_int(gene,nbr_femelles);
        while (femelles[a].patch!=femelles[b].patch || a==b) {
            a = gsl_rng_uniform_int(gene,nbr_femelles);b = gsl_rng_uniform_int(gene,nbr_femelles);
        }
        appa_ff_l = appa_ff_l+apparentement(gene,femelles[a],femelles[b]);

    }
    appa_mm_g = appa_mm_g/nbr_comparaison;
    appa_mf_g = appa_mf_g/nbr_comparaison;
    appa_ff_g = appa_ff_g/nbr_comparaison;
    appa_mm_l = appa_mm_l/nbr_comparaison;
    appa_mf_l = appa_mf_l/nbr_comparaison;
    appa_ff_l = appa_ff_l/nbr_comparaison;

    sortie_apparentement << appa_mm_g <<"\t"<< appa_mf_g <<"\t"<< appa_ff_g <<"\t"<< appa_mm_l <<"\t"<< appa_mf_l <<"\t"<< appa_ff_l <<"\n";
    sortie_apparentement.close();
}
