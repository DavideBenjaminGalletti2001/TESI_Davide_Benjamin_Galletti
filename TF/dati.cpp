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
#include "./LIB/dati.h"

//COSTRUTTORI
Dati::Dati()
{
    //APERTURA FILE INPUT
    fstream input_file {"input.dat", ios::in};
    if (!input_file.good()) throw runtime_error("Errore, non e' possibile aprire il file: input.dat");

    //LETTURA DATI INPUT
    input_file >> pot_chimical_inf;//ESTREMO INFERIORE INTERVALLO DI RICERCA BISEZIONE
    input_file >> pot_chimical_sup;//ESTREMO SUPERIORE INTERVALLO DI RICERCA BISEZIONE
    input_file >> num_points_mesh;//NUMERO PUNTI MESH: INCLUSI I BORDI

    input_file >> L;//LUNGHEZZA INTERVALLO Z
    input_file >> K;//CAMPO ESTERNO
    input_file >> prec;//TOLLERANZA AUTO-CONSISTENZA
    input_file >> acc;//TOLLERANZA BISEZIONE
    input_file >> sigma;//SIGMA
    input_file >> alpha;//MIXING
    n_ave = sigma / L;//DENSITA' MEDIA
    fix_sigma = sigma;//VALORE FISSATO E NON VARIABILE DI SIGMA


    //CREAZIONE MESH
    x_mesh = this->mesh_lin(num_points_mesh, L);  //Punti della mesh totali(logartmica e lineare)

    //NUMERO PUNTI MESH: num_points_mesh - 2
    dim = x_mesh.size();


    //CREO OGGETTO INTEGRATORE
    myT = new trapezi {}; //INTEGRATORE: METODO TRAPEZI
    myBis = new Bisezione {}; //RICERCA ZERI: METODO BISEZIONE
}

void Dati::print_input_dati() //STAMPO A VIDEO DATI INPUT
{
    //STAMPO A VIDEO DATI INPUT
    cout << "%%%%%%%%%%%%%%%%%%" << endl;
    cout << "x_mesh dim = " << dim << endl;
    cout << "pot_chimical_inf = " << pot_chimical_inf << endl;
    cout << "pot_chimical_sup = " << pot_chimical_sup << endl;
    cout << "num_points_mesh = " << num_points_mesh << endl;
    cout << "L = " << L << endl;
    cout << "k = " << K << endl;
    cout << "prec = " << prec << endl;
    cout << "acc = " << acc << endl;
    cout << "sigma = " << sigma << endl;
    cout << "alpha = " << alpha << endl;
    cout << "n_ave = " << n_ave << endl;
    cout << "%%%%%%%%%%%%%%%%%%" << endl;
};
void Dati::print_file_xmesh()//STAMPO SU FILE PUNTI MESH
{
    //MEMORIZZO SU FILE PUNTI MESH
    fstream output_x {"DAT/x.dat", ios::out};
    for (unsigned int i = 0; i < dim; i++) output_x << x_mesh[i] << endl;
}
vector <double> Dati::mesh_lin(unsigned int N, double Max) //CREAZIONE PUNTI MESH LINEARE
{    
    //calcolo min
    double Min = 0;
    
    //calcolo il passo p della mesh lineare
    double p = (Max-Min)/(double)(N-1);
    cout << "Il passo è " << p << endl;
    
    vector<double> vett;
    
    //inserisco mesh lineare
    for(int i=0; i<=N-1; ++i){
    	double val=Min+i*p; 
        vett.push_back(val); 
    }   
    
    return vett;    
}


