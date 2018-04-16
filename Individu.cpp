#ifndef INDIVIDU_CPP_INCLUDED
#define INDIVIDU_CPP_INCLUDED
#include "Individu.h"
#include "patchs.h"
#include <bitset>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>
#include <vector>

#include <iostream>
#include <fstream>

void Individus::initialisation(const gsl_rng* &gene,int nbr_patch,int chromosexe_a,int chromosexe_b,double disp_init_x,double disp_init_y,int taille_bitset){
    etat=true;
    patch=gsl_rng_uniform_int(gene, nbr_patch);
    dispersant=false;
    fecondite = 0;
    progeniture = 0;
    sexe_A=chromosexe_a;
    sexe_B=chromosexe_b;
	if (sexe_A==0) {
		if (disp_init_y > 0) {
			disp_A=disp_init_y+gsl_ran_gaussian(gene,1);
		} else {
			disp_A=-disp_init_y;
		}
	} else if (sexe_A==1) {
		if (disp_init_x > 0) {
			disp_A=disp_init_x+gsl_ran_gaussian(gene,1);
		} else {
			disp_A=-disp_init_x;
		}
	}


	if (sexe_B==0) {
		if (disp_init_y > 0) {
			disp_B=disp_init_y+gsl_ran_gaussian(gene,1);
		} else {
			disp_B=-disp_init_y;
		}
	} else if (sexe_B==1) {
		if (disp_init_x > 0) {
			disp_B=disp_init_x+gsl_ran_gaussian(gene,1);
		} else {
			disp_B=-disp_init_x;
		}
	}

    for (int compteur_loc=0;compteur_loc<taille_bitset;compteur_loc++) {
        neutres_A[compteur_loc]=gsl_rng_uniform_int(gene,2);
        neutres_B[compteur_loc]=gsl_rng_uniform_int(gene,2);
    }
};


void Individus::dispersion(const gsl_rng* &gene,std::vector<Patchs> &parcelles,double &mortalite,int &nbr_patch,int &distance) {
    // moyenne des allèles pour obtenir la stratégie de dispersion
    double disp = 0;
    int vitesse = distance;
    if (distance < 0 ) {
        vitesse = -distance;
    }
	disp = (disp_A+disp_B)/2;

    // test de valeur face à densité du patch
    if ( (distance > 0 && disp < parcelles[patch].densi) || distance < 0) {
        // test de réalisation de la dispersion
        double nbr_alea = gsl_rng_uniform(gene);
        if ( (distance > 0 && (1-(disp/parcelles[patch].densi)) < nbr_alea) || (distance < 0 && nbr_alea < disp) ) {
            // on a un dispersant de plus
            dispersant = true;
            // test de survie du dispersant
            if (mortalite < gsl_rng_uniform(gene)) {
                // le dispersant a survecu
                // attribution d'un nouveau patch au hasard
                if (vitesse>=nbr_patch/2) {
                    int nouveau_patch = gsl_rng_uniform_int(gene, nbr_patch);
                    while (patch==nouveau_patch) {
                        nouveau_patch = gsl_rng_uniform_int(gene, nbr_patch);
                    }
                    patch = nouveau_patch;
                } else { // les individus n'ont pas un déplacmeent illimité
                    // l'environement est circulaire, on peut passer aux patchs
                    // en dessus ou en dessous (droite ou gauche) dans un certain "rayon"
                    // on test si l'individu va en dessou ou au dessus
                    if (gsl_rng_uniform_int(gene, 2)==0) {
                        // dans le negatif
                        int mouvement = -gsl_rng_uniform_int(gene, vitesse)-1;
                        // on test si on descend en dessous de 0 et on ajuste la circularité si nécessaire
                        if (patch+mouvement<0) {
                            patch = patch+mouvement+nbr_patch;
                        } else {
                            patch = patch+mouvement;
                        }
                    } else {
                        // dans le positif
                        int mouvement = gsl_rng_uniform_int(gene, vitesse)+1;
                        // on test si on descend en dessous de 0 et on ajuste la circularité si nécessaire
                        if (patch+mouvement>=nbr_patch) {
                            patch = patch+mouvement-nbr_patch;
                        } else {
                            patch = patch+mouvement;
                        }
                    }
                }
            } else {
                // l'individu est mort pendant sa dispersion
                etat = false;
            }
        } else {
            // l'individu n a pas dispersé, densité haute mais probabilité non réalisée
            dispersant= false;
        }
    } else {
        // l'individu n a pas dispersé car densité trop faible
        dispersant = false;
    }
};

Individus::Individus() {
};


Individus::Individus(const gsl_rng* &gene,Individus &papa, Individus &maman,double &taux_de_mut_s,double &importance_des_mut_s,double &taux_de_mut_n,double &sexratio) {
    // mise à la base des infos non génétiques
    etat = true;
    patch = maman.patch;
    dispersant = false;
    fecondite = 0;
    progeniture = 0;
    // heritabilite du sexe et de la dispersion, via la mere
    // autosomes
    if (gsl_rng_uniform(gene) < sexratio) {
        sexe_A = maman.sexe_A;
		disp_A = maman.disp_A;
    } else {
        sexe_A = maman.sexe_B;
		disp_A = maman.disp_B;
    }


    // mutation (ou pas) de l'allèle de dispersion
    if (gsl_rng_uniform(gene) < taux_de_mut_s) {
        disp_A = disp_A+gsl_ran_gaussian(gene,importance_des_mut_s);
    }

    // autosome
    if (gsl_rng_uniform(gene) < sexratio) {
        sexe_B = papa.sexe_A;
        disp_B = papa.disp_A;
    } else {
        sexe_B = papa.sexe_B;
		disp_B = papa.disp_B;
    }

    // mutation (ou pas) de l'allèle de dispersion
    if (gsl_rng_uniform(gene) < taux_de_mut_s) {
       disp_B = disp_B+gsl_ran_gaussian(gene,importance_des_mut_s);
    }


    // heritabilite des alleles neutres


        neutres_A = maman.neutres_A;
        neutres_B = papa.neutres_A;
    if (taux_de_mut_n >=0) {
        for (int i=0;i<32;i++) {
            if (gsl_rng_uniform_int(gene,2)==0) {
                neutres_A[i] = maman.neutres_B[i];
            }
            if (gsl_rng_uniform_int(gene,2)==0) {
                neutres_B[i] = papa.neutres_B[i];
            }
        }
        // mut neutre A
        if (gsl_rng_uniform(gene) < taux_de_mut_n) {
            // une allele est reversée
            neutres_A.flip(gsl_rng_uniform_int(gene, 32));
        }
        // mut neutre A
        if (gsl_rng_uniform(gene) < taux_de_mut_n) {
            // une allele est reversée
            neutres_B.flip(gsl_rng_uniform_int(gene, 32));
        }
    }


};

#endif // INDIVIDU_CPP_INCLUDED
