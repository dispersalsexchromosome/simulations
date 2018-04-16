#ifndef PATCHS_H_INCLUDED
#define PATCHS_H_INCLUDED
#include <vector>

struct Patchs
{
    double nbr_males;
    double nbr_femelles;
    double densi;
    double sexratio;
    double densi_post_disp;
    double sexratio_post_disp;
    double densijuv;
    double grandlambda;
    double survie;
    std::vector<int> males;
    std::vector<bool> males_use;
    std::vector<int> femelles;

    void ajout_males(int &num_male);
    void ajout_femelles(int &num_femelle);
    void compte_males();
    void retire_male(int &numlocal);
    void use_male(int &numlocal);

};

#endif // PATCHS_H_INCLUDED
