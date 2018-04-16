#include "heterozigotie.h"
#include "Individu.h"
double heterozigotie(Individus &focal) {
    std::bitset <32> resultat = (focal.neutres_A^focal.neutres_B);
    return(double(resultat.count())/32);
}

