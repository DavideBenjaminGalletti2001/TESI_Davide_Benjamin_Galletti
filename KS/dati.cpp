#include "./LIB/dati.h"
Dati::Dati()
{
    //APERTURA FILE INPUT
    cout << "LETTURA FILE INPUT" << endl;
    fstream input_file {"input.dat", ios::in};
    if (!input_file.good()) throw runtime_error("Errore, non e' possibile aprire il file: input.dat");

    //LETTURA DATI INPUT
    input_file >> num_points_mesh;//NUMERO PUNTI MESH: INCLUDE I BOIRD
    input_file >> L;//LATO Z
    input_file >> K;//CAMPO ESTERNO
    input_file >> prec; //TOLLERANZA AUTOCONSISTENZA
    input_file >> sigma; //DENSITA' MEDIA
    input_file >> toll; //TOLLERANZA BISEZIONE
    input_file >> alpha; //MIXING


    
    n_ave =  (sigma/L);//PARTICELLE PER UNITA' DI VOLUME

    //CREAZIONE MESH
    this->mesh_lin(num_points_mesh, L);  //PUNTI DELLA MESH TOTALI(LOGARTMICA E LINEARE)

    //NUMERO PUNTI MESH
    dim = x_mesh.size();//SARA' PARI A num_mesh - 2: CIOE' I BORDI CHE SODDISFANO I VINCOLI DI ANNULLAMENTO

    

}
void Dati::print_input_dati()
{
    //STAMPO A VIDEO DATI INPUT
    cout << "%%%%%%%%%%%%%%%%%%" << endl;
    cout << "x_mesh dim = " << dim << endl;
    cout << "num_points_mesh = " << num_points_mesh << endl;
    cout << "L = " << L << endl;
    cout << "k = " << K << endl;
    cout << "prec = " << prec << endl;
    cout << "n_ave = " << n_ave << endl;
    cout << "toll = " << toll << endl;
    cout << "sigma = " << sigma << endl;
    cout << "%%%%%%%%%%%%%%%%%%" << endl;
};

void Dati::print_file_xmesh()//STAMPA SU FILE PUNTI MESH
{
    //MEMORIZZO SU FILE PUNTI MESH
    fstream output_x {"DAT/x.dat", ios::out};
    for (int i = 0; i < dim; i++) output_x << x_mesh[i] << endl;
}

void Dati::mesh_lin(int N, double Max)
{    
    //calcolo min
    double Min = 0;
    
    //calcolo il passo p della mesh lineare
    diff = (Max-Min)/(double)(N-1);
    cout << "mesh: " << endl;
    cout << "Il passo è " << diff << endl;
    
    //inserisco mesh lineare
    for(int i=1; i<=N-2; ++i){
    	double val=Min +i*diff; 
        x_mesh.push_back(val); 
    } 
    for(int i=0; i<=N-1; ++i){
    	double val=Min +i*diff; 
        x_mesh_boundary.push_back(val); 
    }    
     
}


