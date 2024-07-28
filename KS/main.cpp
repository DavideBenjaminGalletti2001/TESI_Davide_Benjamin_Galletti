#include "LIB/KS.h"


int main()
{
    cout << "BEGIN PROGRAM" << endl;

    KS myKS;

    //STAMPA DATI INPUT E GRANDEZZA DI PARTENZA.
    myKS.print_input_dati();
    myKS.print_file_xmesh();
    myKS.print_file_T();

    //IMPLEMENTO IL CICLO DI KS DI AUTOCONSISTENZA
    myKS.eqKS() ;

    //STAMPA RISULTATI
    //myKS.print_file_eigenfunctions(4);//OUT
    //myKS.print_file_eigenvalues_result();//APP
    myKS.print_file_density();//OUT
    //myKS.print_file_energy();//APP
    //myKS.print_file_subbands();//APP
    //myKS.print_file_sigma_alpha();

    //delete myKS;

    cout << "END PROGRAM" << endl;

    return 0;

}