#include <iostream>
#include <fstream>

#include <vector>
#include <bitset>

#include <cmath>

#include <time.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

#include "Individu.h"
#include "patchs.h"
#include "moyenne.h"
#include "ecart_type.h"
#include "moy_indiv.h"
#include "moy_patch.h"
#include "ecart_type_patch.h"
#include "juv.h"
#include "apparentement.h"
#include "structure_pop.h"
#include "heterozigotie.h"
// main du double expression

int main()
{
    // affichage du lancement du programme
    cout << "init\n";
    //######################################################
    // initialisation des caractéristiques de la simulation
    //######################################################
    //nombre de génération par réplicat //maximal number of generation
    int generation_max=100;
    //nombre de réplicat par jeu de paramètre //number of monte-carlo by parameters set
    int nbrreplicat=2;
    // gestion du sexe heterogametique //who is the heterogametic sex ?
    // sexe 1 = mâle , 0: femelles
    int sexe_heterogametique = 1;
    // ligne qui attribue le sexe_homogametique à l'inverse
    int sexe_homogametique; if (sexe_heterogametique == 1 ) { sexe_homogametique = 0; } else if (sexe_heterogametique==0) { sexe_homogametique = 1;}
    int x=1;int y=0; // rendre plus clair les choses avec du x/y plutôt que des chiffres

    // TYPE DE SIMUL // intensity of the fecundity variation 
    // commenter la solution évitée
    char typesimule[10];
    //sprintf(typesimule,"insecte"); int ts=1;
    sprintf(typesimule,"vertebre"); int ts=2;
    //sprintf(typesimule,"fixe"); int ts=3;
    int sireprofixe=10; //fecundity if fixed

    // NUMERO DE LA SIMUL //simulation number, to change to avoid overwrite
    int numero = 115;
	
    // cout specifique aux femelles //dispersal cost of female
    // si <0, même cout que les males // same as male if <0
    double mu_f = -1;
	
    // valeur d'initialisation spécifique aux homozigotes // allele value for Y
    // si =0, identique à la valeur générales
    double pc_init_e = -0.25;									double pc_e = pc_init_e;
	
     // valeur d'initialisation spécifique aux heterozigotes // allele value for X
    // si =0, identique à la valeur générales
    double pc_init_o = -0.01;									double pc_o = pc_init_o;

	// stochasticité environementale //environmental stochasticity
    unsigned int compteur_sigma=0; double gamme_sigma[] = {1};

    // gestion mélange de la population (casse de la structure génétique) // shuffle
    // 0, aucun melange, 1 melange des deux, 2 melange male, 3 femelles // 0, no shuffle, 1 shuffle, 2, male shuffle, 3, female shuffle
    unsigned int compteur_melange=0; double gamme_melange[] = {0};

    // -2 harem -1 monogamie, 0 monoandrie, 1 polyandrie // mating system
    // probabilité de transition d'un père à un autre pour la stratégie polyandre
    unsigned int compteur_changepere=0; double gamme_changepere[] = {-1};

    // distance de dispertion, si chargé négativement, on a une dispersion à % fixe // dispersal distance, if <0, fixed
    unsigned int compteur_dist_disp=0; int gamme_dist_disp[] = {-200};

    // K //carying capacity
    unsigned int compteur_K ; int gamme_K[] = {100};

    // nbr_patch
    unsigned int compteur_nbr_patch ; int gamme_nbr_patch[] = {100};

    //sex-ratio
    unsigned int compteur_sexratio ; double gamme_sexratio[] = {0.50};

    // lambda //growth rate
    unsigned int compteur_lambda ; double gamme_lambda[] = {2};

    // beta //competition parameter
    unsigned int compteur_beta ; double gamme_beta[] = {1};

    // mu (mortalité à la dispersion) //dispersal mortality/cost
    unsigned int compteur_mu ; double gamme_mu[] = {0.1};

    // in (intensite de l'inbreeding depression, diminution de la fecondite) /inbreeding intensity
    unsigned int compteur_in ; double gamme_in[] = {0};

    // pcinit //si negatif, fixe et pas variant au départ //dispersal initial allele value, everibody have the same if <0, else, randomly drawn around the mean
    unsigned int compteur_pc_init ; double gamme_pc_init[] = {-0.5};

    // taux mut s //mutation rate on dispersal alleles
    unsigned int compteur_taux_muta_s ; double gamme_taux_muta_s[] = {0};

    // force mut s //mutatiion strenght on dispersal alleles
    unsigned int compteur_force_muta_s ; double gamme_force_muta_s[] = {0.05};
	
    // taux mut n, si chargé négativement, les neutres sont annulés //neutral mutation rate
    unsigned int compteur_taux_muta_n ; double gamme_taux_muta_n[] = {0.001};

    //######################################################
    int nbrs_neutres=32;
    // initialisation de trucs pour le suivi
    int taille_simul = sizeof(gamme_sexratio)/sizeof(double)*sizeof(gamme_sigma)/sizeof(double)*sizeof(gamme_changepere)/sizeof(double)*sizeof(gamme_melange)/sizeof(double)*sizeof(gamme_dist_disp)/sizeof(int)*sizeof(gamme_K)/sizeof(int)*sizeof(gamme_nbr_patch)/sizeof(int)*sizeof(gamme_lambda)/sizeof(double)*sizeof(gamme_beta)/sizeof(double)*sizeof(gamme_mu)/sizeof(double)*sizeof(gamme_in)/sizeof(double)*sizeof(gamme_pc_init)/sizeof(double)*sizeof(gamme_taux_muta_s)/sizeof(double)*sizeof(gamme_force_muta_s)/sizeof(double)*sizeof(gamme_taux_muta_n)/sizeof(double)*nbrreplicat;
    int compteur_simul=0;
    int compteur_params=0;
    //######################################################
    // ouverture des fichiers de résultats et de paramètres
    //######################################################


	char nom_parametres[100];
    sprintf(nom_parametres,"parametres_%d.txt",numero);
    std::ofstream parametres(nom_parametres, std::ios::out);

    parametres << "numero"<<"\t"<<"sex_hetero"<<"\t"<<"type_simul"<<"\t"<<"gene"<<"\t"<<"nbrreplicat"<< "\t" <<"si_repro_fixe"<< "\t" << "taille_simul"<< "\t" <<"mu_f"<< "\t" <<"pc_init_o"<< "\t" <<"pc_init_e"<<"\n";
	parametres << numero<<"\t"<<sexe_heterogametique<<"\t"<< typesimule<<"\t"<<generation_max<< "\t"<<nbrreplicat<<"\t"<<sireprofixe<<"\t" <<taille_simul<< "\t" << mu_f << "\t" << pc_init_o << "\t" << pc_init_e <<"\n";
	parametres.close();

    char nom_suivi[100];
    sprintf(nom_suivi,"suivi_%d.txt",numero);
    std::ofstream suivi(nom_suivi,std::ios::out);
    suivi.close();

    //######################################################
    // déclaration des compteurs génétiques
    // et des valeur test des boucles while
    //######################################################
    int compteur1(0);
    //int compteur2(0);
    //int compteur3(0);
    int compteur4(0);
    int replicat(0);
    int compteur_gene(0);
    int compteur_ecriture(0);
    bool taille_pop_ok=true;

    //######################################################
    // déclaration de différentes variables utiles
    //######################################################

    // pour la reproduction et la création des jeunes
    vector<Juv> Jeunes;
    int fecondite=0; // variable de la fécondité d'une femelle
    double hetero_f=0;
    int j_nbr=0; // variable du nombre de jeune totale produit à une génération
    int pere=0; // variable retenant le numero du père
    int mere=0; // variable retenant le numero de la mere
    vector <bool> deja_pere; // vecteur stockant si le mâle est déjà reproduit
    int num_local_pere=0; // variable pour garder la position dans le patch du mal choisi

    //######################################################
    // déclaration et mise en place du générateur aléatoire
    //######################################################
    // creation du generateur
    const gsl_rng *generateur = gsl_rng_alloc(gsl_rng_taus2);
    // init de la graine du générateur
    gsl_rng_set(generateur, time(0));

    // variables liées aux statistiques
    double moyenne_densite_patch=0; // moyenne de densité des patchs (pour la phase adulte)
    double moyenne_sexratio_patch=0; // moyenne de sex-ratio des patchs (pour la phase adulte)
    double moyenne_densite_juv=0; // moyenne de densité des patchs (pour la phase jeune)
    double occupation=0; // taux d'occupation des patchs
    bool ecriture=false; // écriture ou non des résulatts à cette génératio
    int f_dispersant=0; // nbr de dispersant femelle
    int f_survivant=0; // nbr de survivant à la dispersion femelle
    int m_dispersant=0; // nbr de dispersant male
    int m_survivant=0; // nbr de survivant à la dispersion male
    int m_participant=0; //nbr de males se reproduisant
    int f_participant=0; // nbr de femelles différentes se reproduisant
    vector <int> male_dominant; // en cas de harem, vecteur des males dominants
    vector <double> densi_disp_m; // densité subie par les dispersants males
    vector <double> densi_philo_m; // densité subie par les dispersants femelle
    vector <double> densi_disp_f; // densité subie par les philopatriques males
    vector <double> densi_philo_f; // densité subie par les philopatriques femelles


//####################################################################
    // mesures du temps de calcul
    time_t debut = time(NULL); // prise du temps initial
    time_t timedif; // déclaration de la variable du temps au cours de la simulation
    //cout << "\n ****Début des boucles**** \n";
    // pour les différentes stratégies d'appariement possibles
    for (compteur_changepere=0; compteur_changepere < sizeof(gamme_changepere)/sizeof(double); compteur_changepere++) {
    double changepere = gamme_changepere[compteur_changepere];
    // pour les différents degrés de mélange de la population possibles
    for (compteur_melange=0; compteur_melange<sizeof(gamme_melange)/sizeof(double); compteur_melange++) {
    double melange = gamme_melange[compteur_melange];
    // pour les différentes valeur de stochasticité environementale
    for (compteur_sigma=0; compteur_sigma<sizeof(gamme_sigma)/sizeof(double); compteur_sigma++) {
    double sigma = gamme_sigma[compteur_sigma];
    // pour la distance de dispersion
    for(compteur_dist_disp=0;compteur_dist_disp<sizeof(gamme_dist_disp)/sizeof(int);compteur_dist_disp++) {
    int dist_disp = gamme_dist_disp[compteur_dist_disp];
    // pour le K
    for(compteur_K=0;compteur_K<sizeof(gamme_K)/sizeof(int);compteur_K++) {
    int K= gamme_K[compteur_K];
    // pour le nbr de patch
    for(compteur_nbr_patch=0;compteur_nbr_patch<sizeof(gamme_nbr_patch)/sizeof(int);compteur_nbr_patch++) {
    int nbr_patch = gamme_nbr_patch[compteur_nbr_patch];
    // pour le sexratio
    for(compteur_sexratio=0;compteur_sexratio<sizeof(gamme_sexratio)/sizeof(double);compteur_sexratio++) {
    double sexratio = gamme_sexratio[compteur_sexratio];
    // pour lambda
    for(compteur_lambda=0;compteur_lambda<sizeof(gamme_lambda)/sizeof(double);compteur_lambda++) {
    double lambda = gamme_lambda[compteur_lambda];
    // pour beta
    for(compteur_beta=0;compteur_beta<sizeof(gamme_beta)/sizeof(double);compteur_beta++) {
    double beta = gamme_beta[compteur_beta];
    // pour mu
    for(compteur_mu=0;compteur_mu<sizeof(gamme_mu)/sizeof(double);compteur_mu++) {
    double mu = gamme_mu[compteur_mu];
    if (mu_f < 0) {
        mu_f = mu;
    }
    // pour in
    for(compteur_in=0;compteur_in<sizeof(gamme_in)/sizeof(double);compteur_in++) {
    double in = gamme_in[compteur_in];
	
    // pour pcinit
    for(compteur_pc_init=0;compteur_pc_init<sizeof(gamme_pc_init)/sizeof(double);compteur_pc_init++) {
    double pc_init = gamme_pc_init[compteur_pc_init];
	cout << "\n" << pc_init_o << "\n" << pc_init_e << "\n"<< pc_init << "\n";
	if (pc_init_o == 0) {
		pc_o = pc_init;
	}
	if (pc_init_e == 0) {
		pc_e = pc_init;
	}
	
    // pour taux muta s
    for(compteur_taux_muta_s=0;compteur_taux_muta_s<sizeof(gamme_taux_muta_s)/sizeof(double);compteur_taux_muta_s++) {
    double taux_muta_s = gamme_taux_muta_s[compteur_taux_muta_s];
    // pour force muta s
    for(compteur_force_muta_s=0;compteur_force_muta_s<sizeof(gamme_force_muta_s)/sizeof(double);compteur_force_muta_s++) {
    double force_muta_s = gamme_force_muta_s[compteur_force_muta_s];
    // pour beta
    for(compteur_taux_muta_n=0;compteur_taux_muta_n<sizeof(gamme_taux_muta_n)/sizeof(double);compteur_taux_muta_n++) {
    double taux_muta_n = gamme_taux_muta_n[compteur_taux_muta_n];

    //######################################################
    // declaration des différents vecteurs
    // des caractéristiques de la population
    //######################################################
    // nbr d'individus de chaque sexe initialement
    int nbr_i=gamme_K[compteur_K]*gamme_nbr_patch[compteur_nbr_patch];

    // vecteurs initiaux d'un réplicat
    vector <Individus> Males_i(nbr_i);
    vector <Individus> Femelles_i(nbr_i);
    // vecteur de patchs
    vector<Patchs> Parcelles_i(nbr_patch);
    // vecteur d'usage au cours du réplicat
    // initialisés avec les vetcuers initiaux
    int m_nbr=nbr_i;
    int f_nbr=nbr_i;
    vector<Individus> Males(Males_i);
    vector<Individus> Femelles(Femelles_i);
    // vecteur de patchs
    vector<Patchs> Parcelles(nbr_patch);
    // vecteurs temporaires
    vector<Individus> Males_t;
    vector<Individus> Femelles_t;
    vector<Individus> Jeunes_t;
    Individus jeune_temp;


    // calcul de "a" entrant dans le calcul de la survie
    double petita = (pow(lambda,(1/beta))-1)/K;

    // ouverture et gestion des fichiers d ecriture

    compteur_params++;

	char nom_caracteristiques[100];
    sprintf(nom_caracteristiques,"caracteristiques_%d_%d.txt",numero,compteur_params);
    std::ofstream caracteristiques(nom_caracteristiques, std::ios::out);

	caracteristiques <<"sigma"<< "\t"<< "strategie"<<"\t"<<"melange"<<"\t"<<"dist_disp"<<"\t"<<"K"<<"\t"<<"nbr_patch"<<"\t"<<"sexratio"<<"\t"<<"lambda"<<"\t"<<"beta"<<"\t"<<"mu"<<"\t"<<"inb"<<"\t"<<"pc_init"<<"\t"<<"taux_muta_s"<<"\t"<<"force_muta_s"<<"\t"<<"taux_muta_n"<<"\n";
    caracteristiques <<sigma <<"\t"<<changepere<<"\t"<<melange<<"\t"<<dist_disp<<"\t"<<K<<"\t"<<nbr_patch<<"\t"<<sexratio<<"\t"<<lambda<<"\t"<<beta<<"\t"<<mu<<"\t"<<in<<"\t"<<pc_init<<"\t"<<taux_muta_s<<"\t"<<force_muta_s<<"\t"<<taux_muta_n<<"\n";
    caracteristiques.close();

	char nom_compagnon[100];
    sprintf(nom_compagnon,"compagnon_%d_%d.txt",numero,compteur_params);
    std::ofstream compagnon(nom_compagnon, std::ios::out);

    compagnon << "replicat" << "\t" << "generation" <<"\n";
    compagnon.close();
	
	
    char nom_taille_pop[100];
    sprintf(nom_taille_pop,"taille_pop_%d_%d.txt",numero,compteur_params);
    std::ofstream taille_pop(nom_taille_pop, std::ios::out);

    taille_pop << "f_debut"<<"\t"<<"m_debut"<<"\t"<<"f_postdisp"<<"\t"<<"m_postdisp"<<"\t"<<"f_fin"<<"\t"<<"m_fin"<<"\t"<<"juv_prod"<<"\n";
    taille_pop.close();


    char nom_resultats[100];
    sprintf(nom_resultats,"resultats_%d_%d.txt",numero,compteur_params);
    std::ofstream resultats(nom_resultats, std::ios::out);


    resultats <<"f_al_disp_a_e"<<"\t"<<"f_al_disp_b_e"<<"\t"<<"m_al_disp_a_e"<<"\t"<<"m_al_disp_b_e"<<"\t"<<"f_al_disp_a_o"<<"\t"<<"f_al_disp_b_o"<<"\t"<<"m_al_disp_a_o"<<"\t"<<"m_al_disp_b_o"<<"\t"<<"f_dispersant"<<"\t"<<"f_survivant"<<"\t"<<"m_dispersant"<<"\t"<<"m_survivant"<<"\t"<<"m_participant"<<"\t"<<"f_participant"<<"\n";
    resultats.close();

    char nom_valeurs_d_interet[100];
    sprintf(nom_valeurs_d_interet,"valeurs_d_interet_%d_%d.txt",numero,compteur_params);
    std::ofstream valeurs_d_interet(nom_valeurs_d_interet, std::ios::out);

    valeurs_d_interet <<"densite_m_av"<<"\t"<<"densite_et_av"<<"\t"<<"occupation_av"<<"\t"<<"sexratio_av_m"<<"\t"<<"sexratio_av_et"<<"\t"<< "densite_m_ap"<<"\t"<< "densite_et_ap"<<"\t"<<"occupation_ap"<<"\t"<<"sexratio_ap_m"<<"\t"<<"sexratio_ap_et"<<"\t"<<"densite_m_juv"<<"\t"<<"densite_et_juv"<<"\n";
    valeurs_d_interet.close();

    char nom_sortie_apparentement[100];
    sprintf(nom_sortie_apparentement,"apparentement_%d_%d.txt",numero,compteur_params);
    std::ofstream sortie_apparentement(nom_sortie_apparentement,std::ios::out);

    sortie_apparentement << "mm_g_avd"<< "\t"<< "mf_g_avd"<< "\t"<< "ff_g_avd"<< "\t"<< "mm_l_avd"<< "\t"<< "mf_l_avd"<< "\t"<< "ff_l_avd"<<"\t"<<"mm_g_apd"<< "\t"<< "mf_g_apd"<< "\t"<< "ff_g_apd"<< "\t"<< "mm_l_apd"<< "\t"<< "mf_l_apd"<< "\t"<< "ff_l_apd"<< "\n";
    sortie_apparentement.close();


    char nom_parents_lien[100];
    sprintf(nom_parents_lien,"parents_lien_%d_%d.txt",numero,compteur_params);
    std::ofstream parents_lien(nom_parents_lien,std::ios::out);
    char nom_parents_mere[100];
    sprintf(nom_parents_mere,"parents_mere_%d_%d.txt",numero,compteur_params);
    std::ofstream parents_mere(nom_parents_mere,std::ios::out);
    char nom_parents_pere[100];
    sprintf(nom_parents_pere,"parents_pere_%d_%d.txt",numero,compteur_params);
    std::ofstream parents_pere(nom_parents_pere,std::ios::out);
    char nom_parents_survie[100];
    sprintf(nom_parents_survie,"parents_survie_%d_%d.txt",numero,compteur_params);
    std::ofstream parents_survie(nom_parents_survie,std::ios::out);
    char nom_parents_patch[100];
    sprintf(nom_parents_patch,"parents_patch_%d_%d.txt",numero,compteur_params);
    std::ofstream parents_patch(nom_parents_patch,std::ios::out);

    char nom_alleles_lien[100];
    sprintf(nom_alleles_lien,"alleles_lien_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_lien(nom_alleles_lien,std::ios::out);
    char nom_alleles_sexe[100];
    sprintf(nom_alleles_sexe,"alleles_sexe_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_sexe(nom_alleles_sexe,std::ios::out);
    char nom_alleles_A_e[100];
    sprintf(nom_alleles_A_e,"alleles_A_e_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_A_e(nom_alleles_A_e,std::ios::out);
    char nom_alleles_B_e[100];
    sprintf(nom_alleles_B_e,"alleles_B_e_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_B_e(nom_alleles_B_e,std::ios::out);
    char nom_alleles_A_o[100];
    sprintf(nom_alleles_A_o,"alleles_A_o_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_A_o(nom_alleles_A_o,std::ios::out);
    char nom_alleles_B_o[100];
    sprintf(nom_alleles_B_o,"alleles_B_o_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_B_o(nom_alleles_B_o,std::ios::out);
    char nom_alleles_patch[100];
    sprintf(nom_alleles_patch,"alleles_patch_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_patch(nom_alleles_patch,std::ios::out);
    char nom_alleles_fecondite[100];
    sprintf(nom_alleles_fecondite,"alleles_fecondite_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_fecondite(nom_alleles_fecondite,std::ios::out);
    char nom_alleles_progeniture[100];
    sprintf(nom_alleles_progeniture,"alleles_progeniture_%d_%d.txt",numero,compteur_params);
    std::ofstream alleles_progeniture(nom_alleles_progeniture,std::ios::out);

    parents_lien <<"lien_unique"<<"\n";    parents_mere <<"mere"<<"\n";    parents_pere <<"pere"<<"\n";    parents_survie <<"survie"<<"\n";    parents_patch <<"patch"<<"\n";
    parents_lien.close();   parents_mere.close();    parents_pere.close();    parents_survie.close();    parents_patch.close();

    alleles_lien <<"lien_unique"<<"\n";    alleles_sexe <<"sexe"<<"\n";    alleles_A_e <<"A_e"<<"\n";    alleles_B_e <<"B_e"<<"\n";    alleles_A_o <<"A_o"<<"\n";    alleles_B_o <<"B_o"<<"\n";    alleles_patch <<"patch"<<"\n";   alleles_fecondite <<"fecondite"<<"\n";   alleles_progeniture <<"progeniture"<<"\n";
    alleles_lien.close();    alleles_sexe.close();    alleles_A_e.close();    alleles_B_e.close();    alleles_A_o.close();    alleles_B_o.close();    alleles_patch.close();    alleles_fecondite.close(); alleles_progeniture.close();

    // pour chaque réplicat
    for (replicat=1;replicat<=nbrreplicat;replicat++) {
        compteur_simul++;

        char nom_tracker[100];
        sprintf(nom_tracker,"tracker_%d.txt",numero);
        std::ofstream tracker(nom_tracker,std::ios::out);
        tracker.close();

        /*###############################################################
        phase d'initiation du réplicat: début
        ###############################################################*/
        Males=Males_i;
        Femelles=Femelles_i;
        Parcelles=Parcelles_i;

        // initialisation des vecteurs d'usage pour le replicat
        if (sexe_heterogametique==1) {
            // le sexe heterogametique est male (non choix)
            for (compteur1=0;compteur1<nbr_i;compteur1++) {
                Males[compteur1].initialisation(generateur,nbr_patch,x,y,pc_o,pc_e,nbrs_neutres);
                Femelles[compteur1].initialisation(generateur,nbr_patch,x,x,pc_o,pc_e,nbrs_neutres);
            }
        } else {
            // le sexe heterogametique est femelle (choix)
            for (compteur1=0;compteur1<nbr_i;compteur1++) {
                Males[compteur1].initialisation(generateur,nbr_patch,y,y,pc_o,pc_e,nbrs_neutres);
                Femelles[compteur1].initialisation(generateur,nbr_patch,x,y,pc_o,pc_e,nbrs_neutres);
            }
        }

        /*###############################################################
        phase d'initiation du réplicat: fin
        ###############################################################*/
        cout << "****Réplicat initié**** \n";

        compteur_gene=1; // mise à 1 du compteur des générations
        compteur_ecriture=1; // mise à 1 du compteur d'écriture des résultats
        taille_pop_ok=true; // statut "true" que la taille de la pop permet de faire la simulation
        // calcul des tailles de pop
        m_nbr = Males.size();
        f_nbr = Femelles.size();
        j_nbr = 0;

        // tant qu'on a pas atteind la dérnière génération
        // et tant que la population est loin de l'extinction
        // on va dérouler les générations (incrementation en fin de boucle)
        while ( (compteur_gene <= generation_max) && ( taille_pop_ok == true ) ) {
		    //suivi << compteur_gene << "\n";
            // remise à 0 des caractéristiques des parcelles

            std::ofstream tracker(nom_tracker,std::ios::app);
            tracker << "a \n";

            Parcelles = Parcelles_i;
            //##################################
            // condition d'écriture des résultats
            //##################################
            // on regarde si on doit écrire les résultats pour cette génération
            if (compteur_gene==compteur_ecriture) {
                ecriture = true;
                // si oui, on fixe la prochaine génération pour laquelle il faudra écrire
                if (compteur_gene<(generation_max-ceil(generation_max/10))) {
                    //compteur_ecriture = compteur_ecriture+ceil(generation_max/250);
                    compteur_ecriture++;
                } else if (compteur_gene<(generation_max-ceil(generation_max/100))) {
                    //compteur_ecriture = compteur_ecriture+ceil(generation_max/1000);
					compteur_ecriture++;
                } else {
                 compteur_ecriture++;
                }
            } else {
                ecriture = false;
            }
            // forçage ou annulation de l ecriture
            //ecriture = true;
            //##################################
            // ecriture d'une partie des résultats
            //##################################
            if (ecriture==true) {
                timedif = difftime(time(NULL),debut);
                cout << compteur_simul <<" sur "<< taille_simul << "; r: " << replicat << " ; g:"<<  compteur_gene << "; f_nbr: " << f_nbr << " ; m_nbr:  " << m_nbr << " ; j_nbr " << j_nbr<< "\ntemps: " << timedif << "\n" ;
                //cout << "alpha \t";
				
				std::ofstream compagnon(nom_compagnon,std::ios::app);
				compagnon << replicat << "\t" << compteur_gene << "\n";
				compagnon.close();
                std::ofstream resultats(nom_resultats,std::ios::app);
                resultats <<moy_indiv_disp_A_e(Femelles)<<"\t"<<moy_indiv_disp_B_e(Femelles)<<"\t"<<moy_indiv_disp_A_e(Males)<<"\t"<<moy_indiv_disp_B_e(Males)<<"\t"<<moy_indiv_disp_A_o(Femelles)<<"\t"<<moy_indiv_disp_B_o(Femelles)<<"\t"<<moy_indiv_disp_A_o(Males)<<"\t"<<moy_indiv_disp_B_o(Males)<<"\t";
                resultats.close();
                std::ofstream taille_pop(nom_taille_pop,std::ios::app);
                taille_pop << f_nbr<<"\t"<<m_nbr<<"\t";
                taille_pop.close();
            }
            tracker << "b \n";
            //##################################
            // calcul de densite des patchs
            //##################################
            // mise à 0 de la densité pour chaque patch
            for (compteur1=0;compteur1<nbr_patch;compteur1++) {
                Parcelles[compteur1].densi = 0;
            }
            // incrementation de la densité du patch pour chaque male
            for (compteur1=0;compteur1<m_nbr;compteur1++) {
                Parcelles[Males[compteur1].patch].nbr_males++;
            }
            // incrementation de la densité du patch pour chaque femelle
            for (compteur1=0;compteur1<f_nbr;compteur1++) {
                Parcelles[Femelles[compteur1].patch].nbr_femelles++;
            }
            // division par K de la densité
            // et compte du taux d'occupation
            occupation=0;
            for (compteur1=0;compteur1<nbr_patch;compteur1++) {
                Parcelles[compteur1].densi = (Parcelles[compteur1].nbr_males+Parcelles[compteur1].nbr_femelles)/K;
                if (Parcelles[compteur1].densi !=0) {
				// on peut mesurer la sex-ratio
                    Parcelles[compteur1].sexratio = Parcelles[compteur1].nbr_males/(Parcelles[compteur1].nbr_males+Parcelles[compteur1].nbr_femelles);
                 // si il y a au moins un individu sur le patch
                 // on le compte comme occupé
                 occupation++;
                 //cout << Parcelles[compteur1].densi<< "\t";
                } else {
					 Parcelles[compteur1].sexratio = -1;
				}
            }
            tracker << "c \n";
            // écriture des ces valeurs d'interêt des patchs AVANT dispersion (et éventuel mélange)
            if (ecriture==true) {
                moyenne_densite_patch = moy_patch_densi(Parcelles);
                moyenne_sexratio_patch = moy_patch_sexratio(Parcelles);
                std::ofstream valeurs_d_interet(nom_valeurs_d_interet,std::ios::app);
                valeurs_d_interet<<moyenne_densite_patch<<"\t"<< ecart_type_patch_densi(Parcelles,moyenne_densite_patch)<<"\t"<< occupation/nbr_patch<<"\t"<<moyenne_sexratio_patch<<"\t"<<ecart_type_patch_sexratio(Parcelles,moyenne_sexratio_patch) <<"\t";
                valeurs_d_interet.close();
             }
            tracker << "d \n";
			//suivi << "a \n";

            //##################################
            // module de mélange de la pop
            // casse ou non de la structure génétique
            //##################################
            // test pour savoir si l'on mélange la population
            if (melange==0) {
            } else if (melange==1) {
                // si oui
                // mélange des males, convertion de vecteur en tableau pour utiliser la fonction "shuffle"
                int m_patch_mel[m_nbr];
                for (compteur1=0;compteur1<m_nbr;compteur1++) {
                    m_patch_mel[compteur1] = Males[compteur1].patch;
                }
                gsl_ran_shuffle(generateur,m_patch_mel,m_nbr,sizeof(int));
                // reconvertion tableau vers vecteur
                for (compteur1=0;compteur1<m_nbr;compteur1++) {
                    Males[compteur1].patch = m_patch_mel[compteur1];
                }
                // mélange des femelles, convertion de vecteur à tableau pour utiliser la fonction "shuffl
                int f_patch_mel[f_nbr];
                for (compteur1=0;compteur1<f_nbr;compteur1++) {
                    f_patch_mel[compteur1] = Femelles[compteur1].patch;
                }
                gsl_ran_shuffle(generateur,f_patch_mel,f_nbr,sizeof(int));
                 // reconvertion tableau vers vecteur
                for (compteur1=0;compteur1<f_nbr;compteur1++) {
                    Femelles[compteur1].patch = f_patch_mel[compteur1];
                }
            } else if (melange==2) {
                // mélange des males, convertion de vecteur en tableau pour utiliser la fonction "shuffle"
                int m_patch_mel[m_nbr];
                for (compteur1=0;compteur1<m_nbr;compteur1++) {
                    m_patch_mel[compteur1] = Males[compteur1].patch;
                }
                gsl_ran_shuffle(generateur,m_patch_mel,m_nbr,sizeof(int));
                // reconvertion tableau vers vecteur
                for (compteur1=0;compteur1<m_nbr;compteur1++) {
                    Males[compteur1].patch = m_patch_mel[compteur1];
                }
            } else if (melange==3) {
                // mélange des femelles, convertion de vecteur à tableau pour utiliser la fonction "shuffl
                int f_patch_mel[f_nbr];
                for (compteur1=0;compteur1<f_nbr;compteur1++) {
                    f_patch_mel[compteur1] = Femelles[compteur1].patch;
                }
                gsl_ran_shuffle(generateur,f_patch_mel,f_nbr,sizeof(int));
                 // reconvertion tableau vers vecteur
                for (compteur1=0;compteur1<f_nbr;compteur1++) {
                    Femelles[compteur1].patch = f_patch_mel[compteur1];
                }
            }
            tracker << "e \n";
            //##################################
            // module de mélange de la pop: fin
            //##################################

            if (ecriture==true) {
                structure_pop_avd(generateur,Males,Femelles,nom_sortie_apparentement);
            }
            tracker << "f \n";
            //##################################
            // module de la dispersion en elle même
            //##################################
            // test de dispersion -> test de mortalité -> tirage du patch de dispersion
            // pour les males
            m_dispersant=0; //mise à 0 des variables de résultat
            m_survivant=0; //mise à 0 des variables de résultat
            // pour chaque mâle
            for (compteur1=0;compteur1<m_nbr;compteur1++) {
                Males[compteur1].dispersion(generateur,Parcelles,mu,nbr_patch,dist_disp);
                if (Males[compteur1].dispersant==true) {
                    m_dispersant++;
                    if (Males[compteur1].etat==true) {
                        m_survivant++;
                    }
                }
            }
            // pour les femelles
            f_dispersant=0;
            f_survivant=0;
            for (compteur1=0;compteur1<f_nbr;compteur1++) {
                Femelles[compteur1].dispersion(generateur,Parcelles,mu_f,nbr_patch,dist_disp);
                if (Femelles[compteur1].dispersant==true) {
                    f_dispersant++;
                    if (Femelles[compteur1].etat==true) {
                        f_survivant++;
                    }
                }
            }
            tracker << "g \n";
            // retrait des morts
            for (compteur1=0;compteur1<m_nbr;compteur1++) {
                if (Males[compteur1].etat==true) {
                    Males_t.push_back(Males[compteur1]);
                }
            }
            Males = Males_t;
            Males_t.clear();

            for (compteur1=0;compteur1<f_nbr;compteur1++) {
                if (Femelles[compteur1].etat==true) {
                    Femelles_t.push_back(Femelles[compteur1]);
                }
            }
            Femelles= Femelles_t;
            Femelles_t.clear();

            // calcul des tailles de pop
            m_nbr = Males.size();
            f_nbr = Femelles.size();
            tracker << "h \n";
            //##################################
            // fin du module de dispersion
            //##################################
            // si on doit écrire les valeurs de densité etc...APRES la dispersion
            if (ecriture==true) {
                structure_pop_apd(generateur,Males,Femelles,nom_sortie_apparentement);
			}
			// on compte la nouvelle densité d'adultes dans les patchs
			// identique au calcul de densité d'avant dispersion
			// mise à 0 de la densité pour chaque patch
			for (compteur1=0;compteur1<nbr_patch;compteur1++) {
				Parcelles[compteur1].densi_post_disp = 0;
				Parcelles[compteur1].nbr_males = 0;
				Parcelles[compteur1].nbr_femelles = 0;
			}
			// incrementation de la densité du patch pour chaque male
			for (compteur1=0;compteur1<m_nbr;compteur1++) {
				Parcelles[Males[compteur1].patch].nbr_males++;
				//Patchs[m_patch[compteur1]].males(push.back(m_patch)
			}
			// incrementation de la densité du patch pour chaque femelle
			for (compteur1=0;compteur1<f_nbr;compteur1++) {
				Parcelles[Femelles[compteur1].patch].nbr_femelles++;
			}
			// division par K de la densité
			// et compte du taux d'occupation
            occupation=0;
            for (compteur1=0;compteur1<nbr_patch;compteur1++) {
                Parcelles[compteur1].densi_post_disp = (Parcelles[compteur1].nbr_males+Parcelles[compteur1].nbr_femelles)/K;
                if (Parcelles[compteur1].densi_post_disp != 0) {
				// on peut mesurer la sexe-ratio
                    Parcelles[compteur1].sexratio = Parcelles[compteur1].nbr_males/(Parcelles[compteur1].nbr_males+Parcelles[compteur1].nbr_femelles);
                 // si il y a au moins un individu sur le patch
                 // on le compte comme occupé
                 occupation++;
                 //cout << Parcelles[compteur1].densi<< "\t";
                } else {
					Parcelles[compteur1].sexratio = -1;
				}
			}

            tracker << "i \n";
			if (ecriture==true) {

                // écriture des ces valeurs d'interêt des patchs après la dispersion
                moyenne_densite_patch = moy_patch_densi_post_disp(Parcelles);
                moyenne_sexratio_patch = moy_patch_sexratio(Parcelles);
                std::ofstream valeurs_d_interet(nom_valeurs_d_interet,std::ios::app);
                valeurs_d_interet << moyenne_densite_patch << "\t" << ecart_type_patch_densi_post_disp(Parcelles,moyenne_densite_patch) << "\t"<< occupation/nbr_patch<< "\t"<<moyenne_sexratio_patch<< "\t"<<ecart_type_patch_sexratio(Parcelles,moyenne_sexratio_patch)<< "\t";
                valeurs_d_interet.close();

                std::ofstream taille_pop(nom_taille_pop,std::ios::app);
                taille_pop <<f_nbr<<"\t"<<m_nbr<<"\t";
                taille_pop.close();

                std::ofstream resultats(nom_resultats,std::ios::app);
                resultats<<f_dispersant<<"\t"<<f_survivant<<"\t"<<m_dispersant<<"\t"<<m_survivant<<"\t";
                resultats.close();

            }
            tracker << "j \n";
            /*###############################################################
            ###############################################################
            phase de dispersion : fin
            phase de naissance et d'héritabilité : début
            #################################################################
            ###############################################################*/
			//suivi << "b \n";

			// stochasticite environementale pour cette generation
            // on tire au hasard le succès moyen de reproduction sur chaque patch
            for (compteur1=0;compteur1<nbr_patch;compteur1++) {
                if (ts==1) {
                    // Lognormal poethke et al
                    Parcelles[compteur1].grandlambda = gsl_ran_lognormal(generateur,lambda,sigma);
                } else if (ts==2) {
                    // normal revue
                    Parcelles[compteur1].grandlambda = exp(lambda)+gsl_ran_gaussian(generateur,sigma);
                } else if (ts==3) {
                    // sans variation
                    Parcelles[compteur1].grandlambda = 0;
                }
                if (Parcelles[compteur1].grandlambda<0) {
                    Parcelles[compteur1].grandlambda = 0;
                }
            }
            tracker << "k \n";
            //##################################
            // module donnant la fécondité des femelles
            // et attribuant les pères
            //##################################
            m_participant=0; // remise à 0 du nbr de male reproduits
            f_participant=0; // remise à 0 du nbr de femelles reproduites
            // pour chaque patch
            // on établit la liste des pères possibles
            // c'est à dire ceux étant du patch sur lequel on travail
            for (compteur1=0;compteur1<m_nbr;compteur1++) {
                Parcelles[Males[compteur1].patch].ajout_males(compteur1);
            }
            male_dominant.clear();
            // on compte le nbr de males de chaque patch
            for (compteur1=0;compteur1<nbr_patch;compteur1++) {
                Parcelles[compteur1].compte_males();
                //cout << Parcelles[compteur1].nbr_males << endl;
                if (changepere==-2) {
                    //cout << "P: "<<compteur1<<" N: " <<Parcelles[compteur1].nbr_males << "\n";
                    if (Parcelles[compteur1].nbr_males > 0 ) {
                        male_dominant.push_back(Parcelles[compteur1].males[gsl_rng_uniform_int(generateur,Parcelles[compteur1].nbr_males)]);
                        m_participant++;
                    }
                }
            }
            tracker << "l \n";
            // on parcour la liste des femelles et on traite leur repro
            // on randomise la liste des femelles
            int liste_mere[f_nbr];
            for (compteur1=0;compteur1<f_nbr;compteur1++) {
                liste_mere[compteur1]=compteur1;
            }
            gsl_ran_shuffle(generateur,liste_mere,f_nbr,sizeof(int));

            tracker << "m \n";

            for (compteur1=0;compteur1<f_nbr;compteur1++) {
                mere = liste_mere[compteur1];
                int situ_patch = Femelles[mere].patch;
                if (Parcelles[situ_patch].nbr_males > 0 ) {
                    hetero_f = heterozigotie(Femelles[mere]);
                    if (hetero_f < 0.5) {
                       // poisson suivant poethke et al, réduite par l'effet de l'homozigotie
                       fecondite = gsl_ran_poisson(generateur,Parcelles[situ_patch].grandlambda*pow(hetero_f*2,in));
                       if (ts==3) {
                           fecondite = sireprofixe;
                       }
                    } else {
                        // poisson suivant poethke et al, sans reduction de l'homozigotie
                        fecondite = gsl_ran_poisson(generateur,Parcelles[situ_patch].grandlambda);
                        if (ts==3) {
                           fecondite = sireprofixe;
                       }
                    }
                    f_participant++; // la femelle s'est reproduit
                    // on met à jour la densite en juvénile du patch
                    Parcelles[situ_patch].densijuv = Parcelles[situ_patch].densijuv+fecondite;
                    // on attribut un père tiré au hasard au premier jeune de la porté, ou on prend le pere dominant
                    if (changepere != -2) { // pere tire au hasard
                        num_local_pere = gsl_rng_uniform_int(generateur,Parcelles[situ_patch].nbr_males);
                        pere = Parcelles[situ_patch].males[num_local_pere];
                    } else { // pere est le male dominant en mode harem
                        pere = male_dominant[situ_patch];
                    }
                    if (Parcelles[situ_patch].males_use[num_local_pere]==true) {
                        m_participant++; //il a participe a la repro
                        Parcelles[situ_patch].use_male(num_local_pere);  //il ne comptera plus ulterieurement
                        if (changepere == -1) {
                            Parcelles[situ_patch].retire_male(num_local_pere);
                        }
                    }
                    Femelles[mere].fecondite = fecondite;
                    // si la fecondite est supéieur à 0
                    if (fecondite > 0) {
                        // on attribut un père tiré au hasard au premier jeune de la porté
                        Jeunes.push_back(Juv(mere,pere,situ_patch));
                        Males[pere].fecondite++;
                        // si cette fécondité est supérieur à 1 et qu'elle produit donc au moins deux jeunes
                        if (fecondite > 1) {
                            // pour le 2ème jeune jusqu'au dernier
                            for (compteur4=2;compteur4<=fecondite;compteur4++) {
                                // on fait le test pour voir si le père change
                                if (gsl_rng_uniform(generateur) < changepere) {
                                    // si le père change, on tir un nouveau père au hasard
                                    num_local_pere = gsl_rng_uniform_int(generateur,Parcelles[situ_patch].nbr_males);
                                    pere = Parcelles[situ_patch].males[num_local_pere];
                                    if (Parcelles[situ_patch].males_use[num_local_pere]==true) {
                                        m_participant++; //il a participe a la repro
                                        Parcelles[situ_patch].use_male(num_local_pere);  //il ne comptera plus ulterieurement
                                    }
                                }
                                Jeunes.push_back(Juv(mere,pere,situ_patch));
                                Males[pere].fecondite++;
                            }
                        }
                    }
                }
            }
            tracker << "n \n";
            // fin de l'attribution des parents
            // on calcule la survie
            for (compteur1=0;compteur1<nbr_patch;compteur1++) {
                Parcelles[compteur1].survie = 1/pow((1+(petita*(Parcelles[compteur1].densijuv))),beta);
            }
            tracker << "o \n";
            // on retient le nombre de jeunes produit avec la taille du vecteur des numeros de mère
            j_nbr = Jeunes.size();
            // écritures de quelques résultats

            if (ecriture==true) {
                moyenne_densite_juv = moy_patch_densijuv(Parcelles);
                std::ofstream valeurs_d_interet(nom_valeurs_d_interet,std::ios::app);
                valeurs_d_interet << moyenne_densite_juv <<"\t"<< ecart_type_patch_densijuv(Parcelles,moyenne_densite_juv)<<"\n";
                valeurs_d_interet.close();

                std::ofstream resultats(nom_resultats,std::ios::app);
                resultats << m_participant<<"\t"<< f_participant<<"\n";
                resultats.close();
            }
            tracker << "p \n";
			//suivi << "c \n";

            //##################################
            // fin du module de fecondité et de patternité
            // début du module d'hérédité
            //##################################
            for (compteur1=0;compteur1<(j_nbr);compteur1++) {
                if (gsl_rng_uniform(generateur)< Parcelles[Jeunes[compteur1].patch].survie) {
                    Jeunes[compteur1].survie=true;
                    jeune_temp = Individus(generateur,Males[Jeunes[compteur1].pere],Femelles[Jeunes[compteur1].mere],taux_muta_s,force_muta_s,taux_muta_n,sexratio);
                    // incrementation de la progeniture des parents
                    Males[Jeunes[compteur1].pere].progeniture++;
                    Femelles[Jeunes[compteur1].mere].progeniture++;
                    // peu varier selon le sexe heterogametique ou homogametique
                    if (sexe_heterogametique==1) { // le sexe heterogametique sont les males
                        if ((jeune_temp.sexe_A!=jeune_temp.sexe_B)) {
                            Males.push_back(jeune_temp);
                        } else {
                            Femelles.push_back(jeune_temp);
                        }
                    } else {
                        if ((jeune_temp.sexe_A==jeune_temp.sexe_B)) {
                            Males.push_back(jeune_temp);
                        } else {
                            Femelles.push_back(jeune_temp);
                        }
                    }
                }
            }
            tracker << "q \n";
            // ecritures des donnes de comptes de parents et jeunes survivants
            if (ecriture==true && compteur_gene==(generation_max-1)) {

                std::ofstream parents_lien(nom_parents_lien,std::ios::app);std::ofstream parents_mere(nom_parents_mere,std::ios::app);std::ofstream parents_pere(nom_parents_pere,std::ios::app);std::ofstream parents_survie(nom_parents_survie,std::ios::app);std::ofstream parents_patch(nom_parents_patch,std::ios::app);
                for (compteur1=0;compteur1<j_nbr;compteur1++) {
                    // on note les parents de tout les juvs produits et leur survie
                    parents_lien << replicat<< "\n";
                    parents_mere << Jeunes[compteur1].mere<< "\n";
                    parents_pere << Jeunes[compteur1].pere<< "\n";
                    parents_survie << Jeunes[compteur1].survie<< "\n";
                    parents_patch << Jeunes[compteur1].patch<< "\n";
                }
                parents_lien.close();parents_mere.close();parents_pere.close();parents_survie.close();parents_patch.close();

                std::ofstream alleles_lien(nom_alleles_lien,std::ios::app);std::ofstream alleles_sexe(nom_alleles_sexe,std::ios::app);std::ofstream alleles_A_e(nom_alleles_A_e,std::ios::app);std::ofstream alleles_B_e(nom_alleles_B_e,std::ios::app);std::ofstream alleles_A_o(nom_alleles_A_o,std::ios::app);std::ofstream alleles_B_o(nom_alleles_B_o,std::ios::app);std::ofstream alleles_patch(nom_alleles_patch,std::ios::app);std::ofstream alleles_fecondite(nom_alleles_fecondite,std::ios::app);std::ofstream alleles_progeniture(nom_alleles_progeniture,std::ios::app);
                for (compteur1=0;compteur1<m_nbr;compteur1++) {
                    alleles_lien <<replicat<<"\n";
                    alleles_sexe <<"M"<<"\n";
                    alleles_A_e <<Males[compteur1].disp_A_e<<"\n";
                    alleles_B_e <<Males[compteur1].disp_B_e<<"\n";
                    alleles_A_o <<Males[compteur1].disp_A_o<<"\n";
                    alleles_B_o <<Males[compteur1].disp_B_o<<"\n";
                    alleles_patch << Males[compteur1].patch<<"\n";
                    alleles_fecondite << Males[compteur1].fecondite<<"\n";
                    alleles_progeniture << Males[compteur1].progeniture<<"\n";
                }
                for (compteur1=0;compteur1<f_nbr;compteur1++) {
                    alleles_lien <<replicat<<"\n";
                    alleles_sexe <<"F"<<"\n";
                    alleles_A_e <<Femelles[compteur1].disp_A_e<<"\n";
                    alleles_B_e <<Femelles[compteur1].disp_B_e<<"\n";
                    alleles_A_o <<Femelles[compteur1].disp_A_o<<"\n";
                    alleles_B_o <<Femelles[compteur1].disp_B_o<<"\n";
                    alleles_patch << Femelles[compteur1].patch<<"\n";
                    alleles_fecondite << Femelles[compteur1].fecondite<<"\n";
                    alleles_progeniture << Femelles[compteur1].progeniture<<"\n";
                }
                alleles_lien.close();alleles_sexe.close();alleles_A_e.close();alleles_B_e.close();alleles_A_o.close();alleles_B_o.close();alleles_patch.close();alleles_fecondite.close();alleles_progeniture.close();
            }
            tracker << "r \n";
			//suivi << "e \n";
            //mortalité des anciens parents
            // _nbr et f_nbr sont toujours valables
            Males.erase(Males.begin(),Males.begin()+m_nbr);
            Femelles.erase(Femelles.begin(),Femelles.begin()+f_nbr);
            Jeunes.clear();
            Jeunes_t.clear();

            tracker << "s \n";
            m_nbr = Males.size();
            f_nbr = Femelles.size();
            if ((m_nbr+f_nbr) < 100) {
                taille_pop_ok=false;
            }
            compteur_gene++;
            if (ecriture==true) {
                std::ofstream taille_pop(nom_taille_pop,std::ios::app);
                taille_pop << f_nbr<<"\t"<<m_nbr<<"\t"<<j_nbr<<"\n";
                taille_pop.close();
            }
            tracker << "t \n";
            tracker.close();
		//suivi << "f \n";
        }// fin des generation
        Males_t.clear();
        Femelles_t.clear();
        Jeunes.clear();
        Jeunes_t.clear();
		//suivi << " g \n";


        timedif = difftime(time(NULL),debut);
        std::ofstream suivi(nom_suivi,std::ios::app);
        suivi <<compteur_simul<<" sur "<<taille_simul<< " en "<<timedif<< " soit " << double(timedif/3600) <<"heures" <<"\n";
        suivi.close();

    }//fin du replicat
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
    }
	//suivi.close();
    // signal console de fin de simulation
    cout << "return ,  temps: " << difftime(time(NULL),debut);
    return 0;
}
