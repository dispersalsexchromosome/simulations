#include <vector>
#include "patchs.h"

 double moy_patch_densi(std::vector<Patchs> vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    for (int compteur_loc=0; compteur_loc<taille; compteur_loc++) {
        resultat = resultat+vecteur[compteur_loc].densi;
    }
    return(resultat/taille);
}

 double moy_patch_densi_post_disp(std::vector<Patchs> vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    for (int compteur_loc=0; compteur_loc<taille; compteur_loc++) {
        resultat = resultat+vecteur[compteur_loc].densi_post_disp;
    }
    return(resultat/taille);
}

 double moy_patch_densijuv(std::vector<Patchs> vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    for (int compteur_loc=0; compteur_loc<taille; compteur_loc++) {
        resultat = resultat+vecteur[compteur_loc].densijuv;
    }
    return(resultat/taille);
}

 double moy_patch_sexratio(std::vector<Patchs> vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    for (int compteur_loc=0; compteur_loc<taille; compteur_loc++) {
		if (vecteur[compteur_loc].sexratio != -1) {
			resultat = resultat+vecteur[compteur_loc].sexratio;
		} else {
			taille = taille-1;
		}
    }
    return(resultat/taille);
}
