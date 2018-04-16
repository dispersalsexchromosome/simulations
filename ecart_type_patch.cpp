#include <vector>
#include <math.h>
#include "patchs.h"

// fonction ecart type sur trois valeur des patchs
double ecart_type_patch_densi(std::vector<Patchs> vecteur, double moyenne_vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    int compteur=0;
    for (compteur=0; compteur<taille; compteur++) {
        resultat = resultat+sqrt(pow(vecteur[compteur].densi-moyenne_vecteur,2));
    }
return(resultat/taille);
}

double ecart_type_patch_densi_post_disp(std::vector<Patchs> vecteur, double moyenne_vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    int compteur=0;
    for (compteur=0; compteur<taille; compteur++) {
        resultat = resultat+sqrt(pow(vecteur[compteur].densi_post_disp-moyenne_vecteur,2));
    }
    return(resultat/taille);
}

double ecart_type_patch_densijuv(std::vector<Patchs> vecteur, double moyenne_vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    int compteur=0;
    for (compteur=0; compteur<taille; compteur++) {
        resultat = resultat+sqrt(pow(vecteur[compteur].densijuv-moyenne_vecteur,2));
    }
return(resultat/taille);
}

double ecart_type_patch_sexratio(std::vector<Patchs> vecteur, double moyenne_vecteur) {
    double resultat=0;
    int taille = vecteur.size();
    int compteur=0;
    for (compteur=0; compteur<taille; compteur++) {
        resultat = resultat+sqrt(pow(vecteur[compteur].sexratio-moyenne_vecteur,2));

        /*if (vecteur[compteur].nbr_males == 0 && vecteur[compteur].nbr_femelles == 0) {
            taille = taille-1;
        } else {
            resultat = resultat+sqrt(pow(vecteur[compteur].sexratio-moyenne_vecteur,2));
        }*/
    }
return(resultat/taille);
}
