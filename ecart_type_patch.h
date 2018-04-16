#ifndef ECART_TYPE_PATCH_H_INCLUDED
#define ECART_TYPE_PATCH_H_INCLUDED

#include <vector>
double ecart_type_patch_densi(std::vector<Patchs> vecteur, double moyenne_vecteur);

double ecart_type_patch_densi_post_disp(std::vector<Patchs> vecteur, double moyenne_vecteur);

double ecart_type_patch_densijuv(std::vector<Patchs> vecteur, double moyenne_vecteur);

double ecart_type_patch_sexratio(std::vector<Patchs> vecteur, double moyenne_vecteur);

#endif // ECART_TYPE_PATCH_H_INCLUDED
