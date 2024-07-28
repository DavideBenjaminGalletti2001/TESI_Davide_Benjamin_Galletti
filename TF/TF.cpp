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

//COSTRUTTORE
TF::TF() : Dati(), E(dim,0), nNew(dim,0), nOld(dim,0)
{
    
}

//INIZIALIZZAZIONE DENSITA' INIZIALE
void TF::init(void)
{
    for(unsigned int i=0; i<dim; i++) nOld[i] = n_ave;//GENERO VETTORE DI CAMPIONAMENTO DISCRETO DEI VALORI DELLA DENSITA' DI PROVA
}


//STUDIO DIPENDENZA DALLA ENERGIA DI FERMI DI SIGMA. STUDIO PER DEFINIRE INTERVALLO DI INTEGRAZIONE
void TF::check_interval(void)
{
    fstream out;
    out.open("DAT/N(mu)_K=" +to_string(K) + ".dat", ios::app);//APERTURA FILE DI OUTPUT E DENOMICAZIONE
    pot_chimical=pot_chimical_inf;//INIZIALIZZO MU AL VALORE MINIMO
    double passo_mu=abs(pot_chimical_sup-pot_chimical_inf)/20;//DEFINISCO IL PASSO DI INCREMENTO CON CUI AUMENTARE MU
    cout << "Passo potenziale chimico = " << passo_mu << endl;//STAMPO PASS INCREMENTO MU
    do{
        init();//INIZIALIZZO LA DENSITA' DI PROVA
        double a=operator()(pot_chimical);//CALCOLO IL NUMERO DI ELETTRONI A FISSATO MU
        out << pot_chimical << " " << a << endl; //STAMPO MU E CORRISPONDENTE SIGMA
        std::cout << pot_chimical<< " " << a << endl;//STAMPO MU E CORRRISPONDENTE SIGMA 
        
        pot_chimical+=passo_mu;
    } while(pot_chimical <= pot_chimical_sup);
    
    out.close();
}

void TF::set_Ekin()//ENERGIA CINETICA
{
    for(int k = 0; k < dim; k++) E[k] = pow(nOld[k],5./3.);
    E_kin = Ctf*myT->integra(E,x_mesh);
}

void TF::set_Eext()//ENERGIA ESTERNA
{
    for (int k = 0; k < dim; k++) E[k] = x_mesh[k]*nOld[k];
    E_ext = K*myT->integra(E,x_mesh);
}
void TF::set_Ecoul()//ENERGIA COULOMBIANA
{
    vector<double> I (dim, 0.);
    for (int k = 0; k < dim; k++){
        for (int j = 0; j < dim; j++) I[j] = (nOld[j]-n_ave)*abs(x_mesh[j]-x_mesh[k]);
        E[k] = (nOld[k]-n_ave)*myT->integra(I,x_mesh);
        //cout << E[k] << endl;
        
    }
    
    E_coul = -M_PI*myT->integra(E, x_mesh);
}

void TF::set_Ex()//ENERGIA DI SCAMBIO
{
    for (int k = 0; k < dim; k++) E[k] = pow(nOld[k],4./3.);
    E_x = -Ctfx*myT->integra(E,x_mesh);
}

void TF::set_Err()//ERRORE AUTO-CONSISTENZA
{
    for (int k = 0; k < dim; k++) E[k] = pow(nOld[k]-nNew[k],2);
    err = myT->integra(E,x_mesh);
}


//RISOLUZIONE EQUAZIONE DI TF A DATO MU
double TF::operator()(double  mu)//CALCOLA NUMERO PARTICELLE DATO MU RISOLVENDO IN MODO AUTOCONSISTENTE L'EQUAZIONE DI THOMAS-FERMI
{
    cout << endl;
    fstream output_fp {"DAT/fp.dat", ios::out};
    
    
    double E_old;
    double E_new;
    init(); //Ad ogni ciclo devo dire la prima densità di prova.
    int counter=0;
    int iter=0;
    
    do{
        counter++;
          
        iter++;
        //fileout << endl << "# " << iter << endl;
        double sum=0.;
        E_old = E_new; 
     
        ofstream fout {"DAT/density.dat"};
        //fill(nNew.begin(),nNew.end(),0);
        for (r = 0; r < dim; r++){
            vector<double> I(dim,0);
            for (int k = 0; k < dim; k++) I[k] = -2*M_PI*(nOld[k]-n_ave)*abs(x_mesh[k]-x_mesh[r]);
            double b=(mu-K*x_mesh[r]-e*e*myT->integra(I, x_mesh) + 4./3.*Ctfx*pow(nOld[r],1/3));
            //cout << potExt(r) << endl;
            //cout << b << endl;
            
            if(b>0){
                nNew[r]=(pow(3./(5.*Ctf)*b,1.5));
             
            } else {
                nNew[r]=0;
            }
            fout << x_mesh[r] << " " <<  nNew[r] << endl;
   
        }
       

        
        
        
        
        for(int i=0; i<dim; i++){
            nOld[i]=(1-alpha)*nOld[i]+alpha*nNew[i]; //Il mixing va bene anche da quando b<0. Posso cambiare i parametri, influisce sulla velocità del programma.
            //La nuova densità diventa quella vecchia
        }
        
        set_Err();
        sigma =  myT->integra(nOld,x_mesh);
        n_ave = sigma /L;
        if(counter%100 == 0) cerr << iter << ") POT. CHIMICO: " << mu << "; " << " ERRORE: " << err << "; SIGMA: " << sigma << endl; //Commentare quando tutto funziona.
      
        set_Ekin();
        set_Eext();
        set_Ecoul();
        set_Ex();
        E_new = E_kin + E_ext +  E_coul+E_x; //energia con la densità nuova.
       
        E_tot = E_kin + E_ext +  E_coul+E_x;
        //print_energy_display();
    } while(err > prec); //Metti nel file input.
    cout << "Errore = " << abs(E_new-E_old) << endl;
    cout << "Ekin = " << E_kin << endl;
    cout << "Eext = " << E_ext<< endl;
    cout << "Ecoul = " << E_coul << endl;
    cout << "Ex = " << E_x << endl;
    cout << "Etot = " << E_tot << endl;
    return myT->integra(nOld, x_mesh);
}


//CALCOLO DENSITA' A DATA SIGMA

void TF::density(void)
{
    
    
    ostringstream oss;
    oss << setprecision(2)<< K; 
    string a_K = oss.str();
    oss.str("");
    oss << setprecision(2) << fix_sigma;
    string a_sigma = oss.str();
    oss.str("");


    string title = "DAT/density_K=" +a_K+"_sigma=" + a_sigma + ".dat";
    fstream fileout (title.c_str(), ios::out);
    
    //fileout.open("DAT/density_K=" +to_string(K)+"_sigma=" + to_string(fix_sigma) + ".dat", ios::out);
    //cout << "f(a) =" << bis(pot_chimical_inf) << endl;
    //cout << "f(b) =" << bis(pot_chimical_sup) << endl;
    
    double root = myBis->esegui(pot_chimical_inf, pot_chimical_sup, acc, bis);//IMPLEMENTO IL METODO DI BISEZIONE NELL'INTERVALLO [muInf,muSup]
                                                                //sigma INDICA DENSITA' SUPERFIICIALE DI PARTICELLE DESIDERATO
                                                                //Bis LA FUNZIONE PER RICAVARE N PARTICELLE A DATO POTENZIALE CHIMICO
                                                                //acc INVECE INDICA LA TOLLERANZA SUL POTENZIALE per
                                                                //DIRE DI AVER TROVATO UNO ZERO
    
    //OTTENIAMO E SALVIAMO IN root IL POTENZIALE INDICANTE LO ZERO
    cout << " 	*** RADICE ***	" << root << endl;
    //Memorizzo la densità associata al potenziale chimico trovato.
    //Iniz();
    double sigma_found=operator()(root);//
    cout << "NUMERO SIGMA TROVATO :  " << sigma_found << endl;
    cout << "Stampo il file" << endl;
    for(int i=0; i<dim; i++)
        fileout << x_mesh[i] << " " << nNew[i] <<  endl;
         
    fileout.close();
}

void TF::print_energy_file(){
    
    
    ostringstream oss;
    oss << setprecision(2)<< K; 
    string a_K = oss.str();
    oss.str("");
    string title = "DAT/energy_K="+a_K+".dat";
    fstream fileout (title.c_str(), ios::app);

    fileout << sigma << " " << E_tot << " " << E_kin << " " << E_ext << " " << E_coul << " " << E_x << " " << errore <<  endl;  
}

void TF::print_energy_display()
{
    cout << "ENERGIA" << endl;
    cout << "Ekin     " << E_kin<< endl;
    cout << "PotExt   " << E_ext << endl;
    cout << "Hartree  " << E_coul << endl;
    cout << "scambio  " << E_x << endl;
    cout << "ETot     " << E_tot << endl;
}