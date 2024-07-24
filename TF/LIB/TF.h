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
#ifndef __TF__
#define __TF__
#include "dati.h"
#include <functional>
#include <sstream>

using namespace std;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
class TF : public Dati
{
    public:
        //------//
        //COSTRUTTORE
        TF();

        void set_Ekin();
        void set_Eext();
        void set_Ecoul();
        void set_Ex();
        void set_Err();

        //------//
        //INIZIALIZZO DENSITA' DI INPUT E DISTANZA MESH:
        void init(void);

        //------//
        //STAMPO NUMERO PARTICELLE IN FUNZIONE DEL POTENZIALE. ALL'INTERNO IMPLEMENTO IL METODO AUTOCONSISTENTE PER LA
        //RISOLUZIONE DELLA EQUAZIONE DI THOMAS-FERMI.
        void check_interval(void);

        //------//
        //METODO PER CALCOLARE IL NUMERO DI PARTICELLE DATO IL POTENZIALE DI INPUT. RISOLUZIONE AUTOCONSISTENTE EQUAZIONE THOMAS-FERMI
        double operator()(double); 

        //------//
        //METODO BISEZIONE PER CALCOLARE LA DENSITA' A DATO NUMERO DI ELETTRONI
        void density(void);

        //------//
        //METODO BISEZIONE PER CALCOLARE LA DENSITA' A DATO NUMERO DI ELETTRONI
        void print_energy_file(void);

        //------//
        //METODO BISEZIONE PER CALCOLARE LA DENSITA' A DATO NUMERO DI ELETTRONI
        void print_energy_display(void);


        
        //9)NUMERO PARTICELLE FOUNDED - NUMERO PARTICELLE DESIDERATO
        std::function<double (double)> bis = [&] (double mu) -> double{
             return operator()(mu) -fix_sigma;
        };

        vector<double> E;
        vector<double> nOld;
        vector<double> nNew;


       
        
        
};
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
#endif