#ifndef __KS__
#define __KS__

//LIBRERIE PERSONALI
#include "dati.h"


//LIBERIE DI SISTEMA
#include <vector>
#include <fstream>
#include <functional>
#include <iomanip>
#include <string>
#include <iomanip>
#include <sstream>

//USING
using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//

class KS : public Dati{
    public:
        //COSTRUTTORE
        KS(); 

        //DISTRUTTORE
        ~KS(){delete [] KS_matrix; delete [] eigenvalues; delete [] work; delete [] ipiv;};

        //PRINT MATRICES ON FILES
        void print_file_T();

        //PRINT OBSERVABLES ON FILES
        void print_file_eigenvalues_result();
        void print_file_density();
        void print_file_energy();
        void print_file_subbands();
        void print_file_eigenvalues();
        void print_file_KS_matrix(int);
        void print_file_eigenfunctions(int);
        void print_file_sigma_alpha();
 
        void eqKS();//AUTOCONSISTENZA
           
    private:
        
        //INIZIALIZZAZIONE MATRICI

        void init_T(); //MATRICE CINETICA
        void set_KS_matrix(); //HAMILTONIANA DI KS
        void set_KS_potential(); //POTENZIALE


        //Restituisce la posizione dell'elemento di matrice (i,j) all'interno dell'array. Nel mio codice le matrici sono trattate come semplici array per avere compatibilità con LAPACK. N.B: i e j da 0 a (N-1)
        int ind(int i, int j) {
            return i*dim+j; //riga i e colonna j. N va ad indicare la lunghezza di una riga
        }

        void density();//IMPLEMENTO NUOVA DENSITA'
        void norm();//NORMALIZZO AUTOFUNZIONI
        void set_init_density();//SETTO DENSITA' INIZIALE
        void mixing(vector<double>& vNew, vector<double>& vOld);//MIXO VECCHIO E NUOVA VALUTAZIONE OSSERVABILE 

        double electron_counting(double ef);//CONTEGGIO ELETTRONI A DATA ENERGIA DI FERMI

        void set_exc();//POTENZIALE SCAMBIO E CORRELAZIONE
        void set_Vxc();//POTENZIALE DERIVATA SCAMBIO E CORRELAZIONE
        void set_Vh();//POTENZIALE HARTREE
        void set_Vext();//POTENZIALE ESTERNO
        


        //%%%%%%%%ENERGIE%%%%%%%%%//

        void set_Eh();//ENERGIA DI HARTRRE
        void set_EH();
        void set_Ecoul();
        void set_Exc();//ENERGIA DI SCAMBIO-CORRELAZIONE: PARAMETRIZZAZIONE PERDEW-YANG
        void set_EVxc();//ENERGIA DERIVATA SCAMBIO E CORRELAZIONE
        void set_Eext();//ENERGIA ESTERNA
        void set_Et();//ENERGIA LIVELLI
        void set_Ekin();//ENERGIA CINETICA
        void set_Eb();

        //%%%%%%%%ERRORE%%%%%%%%%//

        void set_err();//ERRORE
        


        //%%%%%%%%DENSITA' SUPERFICIALE%%%%%%%%%//
        //FUNZIONE PER TROVARE DENSITA' SUPERFICIALE: SI SFRUTTA IL METODO DI BISEZIONE

        std::function<double (double)> fsigma = [&] (double ef_value) -> double{
            return electron_counting(ef_value) - sigma;
        };

        //ARRAY PER LIBRERIA 
        double*  KS_matrix;
        double* eigenvalues;
        

        //LAPACK
        int lwork_dsyev;
        int lwork_dsysv;
        int lwork;
        double *work ;
        int *ipiv ;//Per dsysv 
        int info;//Per dsysv 

        //trapezi
        trapezi tr;
        //bisezione
        Bisezione bi;

        
        vector<double> T;//MATRICE TRIDIAGONALE ENERGIA CINETICA.
        vector<double> Vnew;//NUOVO POTENZIALE
        vector<double> Vold;//VECCHIO POTENZIALE
        vector<double> nNew;//NUOVA DENSITA' 
        vector<double> nOld;//VECCHIA DENSITA'
        vector<double> Exc;//ENERGIA SCAMBIO-CORRELAZIONE
        vector<double> Vh;//ENERGIA DI HARTREE
        vector<double> Vext;//ENERGIA ESTERNA
        vector<double> Vxc;//POTENZIALE SCAMBIO CORRELAZIONE
        vector<double> E;//ENERGIA
        vector<double> Err;//ERRORE


        
        


};

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//



#endif 