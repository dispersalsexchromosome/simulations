#include <vector>
#include "Individu.h"

 double moy_indiv_disp_A(std::vector<Individus> vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    for (int compteur_loc=0; compteur_loc<taille; compteur_loc++) {
        resultat = resultat+vecteur[compteur_loc].disp_A;
    }
    return(resultat/taille);
}

 double moy_indiv_disp_B(std::vector<Individus> vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    for (int compteur_loc=0; compteur_loc<taille; compteur_loc++) {
        resultat = resultat+vecteur[compteur_loc].disp_B;
    }
    return(resultat/taille);
}
