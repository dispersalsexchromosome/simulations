#include <vector>
#include <math.h>
using namespace std;
// fonction ecart type
double ecart_type(vector<double> vecteur, double moyenne_vecteur) {
     double resultat=0;
     int taille = vecteur.size();
     int compteur=0;
     for (compteur=0; compteur<taille; compteur++) {
        resultat = resultat+sqrt(pow(vecteur[compteur]-moyenne_vecteur,2));
     }
     resultat = resultat/taille;
return(resultat);
}
