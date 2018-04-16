#include <vector>
#include "patchs.h"

void Patchs::ajout_males(int &num_male) {
    males.push_back(num_male);
    males_use.push_back(true);
}

void Patchs::ajout_femelles(int &num_femelle) {
    femelles.push_back(num_femelle);
}

void Patchs::compte_males() {
    nbr_males = males.size();
}

void Patchs::retire_male(int &numlocal) {
    males.erase(males.begin()+numlocal);
    males_use.erase(males_use.begin()+numlocal);
    nbr_males = nbr_males-1;
}

void Patchs::use_male(int &numlocal) {
    males_use[numlocal]=false;
}
