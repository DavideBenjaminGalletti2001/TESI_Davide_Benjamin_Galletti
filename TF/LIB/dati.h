/*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
$$$$$$$$$$$$                 $$$$$$$$$$$$               $$$$$$$$$$$$$$$
$$        $$$$               $$        $$$$           $$$$
$$          $$$$             $$          $$$$       $$$$
$$            $$$$           $$            $$$$     $$
$$              $$$$         $$            $$$$     $$
$$                $$$$       $$          $$$$       $$
$$                  $$$$     $$$$$$$$$$$$$$$$$      $$
$$                  $$$$     $$             $$$$    $$        $$$$$$$$
$$                  $$$$     $$             $$$$    $$              $$$$
$$                $$$$       $$           $$$$      $$$$            $$$$
$$              $$$$         $$         $$$$          $$$$        $$$$
$$$$$$$$$$$$$$$$$            $$$$$$$$$$$$$              $$$$$$$$$$$$ 

AUTORE: Davide Benjamin Galletti
DATA: 29/06/2024
UNIVERSITA': Università degli Studi di Milano
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef __dati__
#define __dati__

#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>

#include "trapezi.h"
#include "bisezione.h"

using namespace std;

class Dati{
    
    public:
        //COSTRUTTORE
        Dati();

        //METODI
        vector<double> mesh_lin(unsigned int , double);//costruzione mesh lineare
        void print_input_dati();
        void print_file_xmesh();

        //DATI INPUT 
        double pot_chimical_inf;//valore min pot. chimico  
        double pot_chimical_sup;//valore max pot. chimico 
        unsigned int num_points_mesh;//numero punti mesh 
        double L;//larghezza intervallo asse x
        double K;//campo esterno
        double prec;//tolleranza errore
        double acc;//accettanza bisezione
        double sigma;//densita' superficiale
        double alpha; //mixing
        double fix_sigma;

        //ALTRE VARIABILI
        double n_ave;//densita' su unita' di volume
        int boolean;
        unsigned int dim;
        double corrFactor;
        double errore;
        double Np;
        unsigned int r;
        double pot_chimical;//potenziale chimico
        double E_tot;
        double E_kin;
        double E_x;
        double E_ext;
        double E_coul;
        double err;

        //COSTANTI
        double const mass_electron = 1;//massa elettrone
        double const reduced_planck_const = 1;//costante di planck tagliata
        double const a0 = 1;//raggio atomico
        double const e = 1;//constante e^2 = 1/(4*M_PI*EPSILON0)
        double const Ctfx=3./4.*pow(3/M_PI,1./3.)*e*e;
        double const Ctf=3./5.*pow(2*M_PI*M_PI,2./3.)*reduced_planck_const*reduced_planck_const/(2.*mass_electron);

        //VECTOR
        vector<double> nOld;//vecchia densita'
        vector<double> nNew;//nuova' densita'
        vector<double> x_mesh;//punti della mesh

        //INTEGRATORE
        trapezi* myT;

        //SOLUTORE
        Bisezione* myBis;

        
};

#endif