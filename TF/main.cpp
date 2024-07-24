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


//LIBRERIE PERSONALI
#include "./LIB/TF.h"

//LIBRERIE DI SISTEMA
#include <cstdlib>
#include <iomanip>
#include <iostream>

//USING NAMESPACE
using namespace std;

int main(int argc , char **argv)
{
    cout << "BEGIN PROGRAM" << endl;
    TF myTF;
    myTF.print_input_dati();
    myTF.print_file_xmesh();
    //myTF.check_interval();
    //cout << myTF(5) << endl;
    myTF.density();
    //myTF.print_energy_display();
    cout << "END PROGRAM" << endl;    
    
    return 0;
}