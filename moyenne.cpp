#include <vector>
using namespace std;
// fonction moyenne
double moyenne(const vector<double> vecteur) {
     double resultat=0;
     int taille = vecteur.size();
     int compteur=0;
     for (compteur=0; compteur<taille; compteur++) {
        resultat = resultat+vecteur[compteur];
     }
     resultat = resultat/taille;
return(resultat);
}
