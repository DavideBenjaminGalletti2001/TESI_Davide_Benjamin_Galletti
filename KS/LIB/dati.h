
#ifndef __dati__
#define __dati__

//LIBRERIE DI SISTEMA
#include <vector>
#include <fstream>
#include <cmath>
#include <iostream>


//LIBRERIE PERSONALI
#include "nr.h"
#include "lapack.h"
#include "trapezi.h"
#include "bisezione.h"



using namespace std;

class Dati{
    
    public:
        //COSTRUTTORE
        Dati();

        
        //~Dati(){delete tr; delete bi; delete work; delete ipiv; delete [] KS_matrix; delete [] eigenvalues;};
        
        

        //METODI
        void mesh_lin(int , double);//costruzione mesh lineare
        void print_input_dati();
        void print_file_xmesh();

        //DATI INPUT 

        int num_points_mesh;//numero punti mesh 
        double L;//larghezza intervallo asse x
        double K;//campo esterno
        double prec;//tolleranza autoconsistenza
        double n_ave;//densita' media
        double toll;//tolleranza bisezione
        double alpha;
        

        //ALTRE VARIABILI
        int dim; //numero punti mesh
        double eF; //energia di fermi
        double nF; //numer quantico di fermi
        double kF;  //numero d'onda di fermi
        double sigma; //densita' superficiale
        int r = 0;
        double m_sigma = 0;
        double m_n_ave = 0;

        double Eold = 0;
        double Enew = 0;

        double diff;
         

       

        //COSTANTI
        double const mass_electron = 1;//massa elettrone
        double const reduced_planck_const = 1;//costante di planck tagliata
        double const a0 = 1;//raggio atomico
        double const e = 1;//constante e^2 = 1/(4*M_PI*EPSILON0)
        double const B = mass_electron/(M_PI*reduced_planck_const*reduced_planck_const);
        double const C = (reduced_planck_const*reduced_planck_const/mass_electron);
        

        vector<double> x_mesh;//punti della mesh
        vector<double> x_mesh_boundary;

        double E_h;//ENERGIA DI HARTRRE
        double E_coul;
        double E_H;
        double E_t = 0.;//SOMMA AUTOVALORI
        double E_xc;//ENERGIA SCAMBIO-CORRELAZIONE
        double E_kin;//ENERGIA CINETICA
        double E_ext;//ENERGIA ESTERNA
        double E_eff;//ENERGIA EFFICACE
        double E_e;
        double E_bb;
        double E_b;
        double E_Vxc;//ENERGIE: HARTREE, SINGOLA PARTICELLA, SCAMBIO CORRELAZIONE, SCAMBIO CORRELAZIONE DERIVATA

        double sigma_alpha1 = 0;
        double sigma_alpha2 = 0;
        double sigma_alpha3 = 0;


        double const A1 = 0.9164/2.;
        double const A2 = 0.2846/2.;
        double const A3 = 1.0529;
        double const A4 = 0.3334;
        double const A5 = 0.0960/2.;
        double const A6 = 0.0622/2.;
        double const A7 = 0.0232/2.;
        double const A8 = 0.0040/2.;


        double cutoff = 0.000000000000000001;

   
        




        
};

#endif