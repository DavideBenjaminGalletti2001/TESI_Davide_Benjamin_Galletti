
#include "LIB/KS.h"


KS::KS(): Dati() , tr(), bi(), T(dim*dim, 0), Vnew(dim,0), Vold(dim,0), 
nNew(num_points_mesh, 0), nOld(num_points_mesh,0), 
E(num_points_mesh,0), 
Err(num_points_mesh), 
Exc(dim,0), Vh(dim,0), Vext(dim,0), Vxc(dim,0)              
{
    
    //SETTAGGIO MATRICE KS
    KS_matrix = new double[dim*dim];//matrice da diagonalizzare: definisce l'hamiltonia della equazione di KS. E' linearizzata
    eigenvalues = new double[dim];
    
    //LAPACK: DIAGONALIZZAZIONE MATRICE
    lwork_dsyev = (NB_dsytrd(dim)+2)*dim;
    lwork_dsysv = (NB_dsytrf(dim))*dim;
    lwork = max(lwork_dsyev,lwork_dsysv);
    work = new double[lwork];
    //Per dsysv
    ipiv = new int[dim];
    
    E_bb = -n_ave*n_ave*e*e*M_PI*pow(L,3)/3;
    init_T();//INIZIALIZZIAMO LA MATRICE CINETICA DIRETTAMENTE NEL COSTRUTTORE: E' FISSA E DEFINITA UNA VOLTA PER TUTTE.
    
}


//INIZIALIZZAZIONE MATRICE CINETICA:
void KS::init_T()
{
    //la funzione ind mi permette di interpretare matrici linearizzate come matrici effettive
    T[ind(0,0)] = -2./((x_mesh[1]-x_mesh[0])*(x_mesh[0]));
    T[ind(dim-1,dim-1)]= -2./((L-x_mesh[dim-1])*(x_mesh[dim-1]-x_mesh[dim-2]));
    
    T[ind(0,1)] = +2./((x_mesh[1]-x_mesh[0])*(x_mesh[1]));
    T[ind(1,0)] = +2./((x_mesh[1]-x_mesh[0])*(x_mesh[2]-x_mesh[0]));
                    
    T[ind(dim-2,dim-1)]= +2./((x_mesh[dim-1]-x_mesh[dim-2])*(x_mesh[dim-1]-x_mesh[dim-3]));
    
    T[ind(dim-1,dim-2)]= +2./((x_mesh[dim-1]-x_mesh[dim-2])*(L-x_mesh[dim-2]));
             
             
    for(int i=1; i<=dim-2; ++i){
        T[ind(i,i)] = -2./((x_mesh[i+1]-x_mesh[i])*(x_mesh[i]-x_mesh[i-1]));
    }
    
    for(int i=1; i<=dim-3; ++i){
        T[ind(i,i+1)] = +2./((x_mesh[i+1]-x_mesh[i])*(x_mesh[i+1]-x_mesh[i-1]));
        T[ind(i+1,i)] = +2./((x_mesh[i+1]-x_mesh[i])*(x_mesh[i+2]-x_mesh[i]));        
    }
}

void KS::set_init_density()
{
    
    for (int i = 1; i < dim+1; i++){
        nOld[i]=n_ave;//setto uguale al valor medio per generalita', ma non e' particolarmente influente
        
    }
    
}

void KS::set_exc()
{
    for (int k = 0; k < dim; k++)
    {
        double rs1 = (pow(3/(4*M_PI*nOld[k+1]*a0*a0*a0),1./3.));
        if (rs1 >= 1)
        {
            Exc[k]= -A1/rs1-A2/(1+A3*pow(rs1,1./2.)+A4*rs1);
        }else{
            Exc[k] = -A1/rs1 -A5+A6*log(rs1)-A7*rs1 + A8*rs1*log(rs1);
        }
         
    }
}


void KS::set_Vh()
{
    for (int k = 0; k < dim; k++)
        Vh[k] = -2*M_PI*e*e*(nOld[k+1]-n_ave)*abs(x_mesh[k] - x_mesh[r]);
};

void KS::set_Vext()
{
    for (int k = 0; k < dim; k++)
        Vext[k] = K*(x_mesh[k]);
};



void KS::set_Vxc()
{
    for (int k = 0; k < dim; k++)
    {      
        double rs1 = (pow(3./(4.*M_PI*nOld[k+1]*a0*a0*a0),1./3.));
        double S = 1/(pow(2.*2.*2.*2.*3.*M_PI*M_PI,1./3.))*1/rs1;
        if (rs1 >= 1)
        {
            Vxc[k]= -A1/rs1-A2/(1+A3*pow(rs1,1./2.)+A4*rs1)-
            +S*(A1/(2.*rs1*rs1)+A2/pow(1+A3*pow(rs1,1./2.)+A4*rs1,2.)*(+A3/2.*pow(rs1,-1./2.)+A4));
        }else{
            Vxc[k] = -A1/rs1 -A5+A6*log(rs1)-A7*rs1 + A8*rs1*log(rs1)-
            +S*(A1/(2.*rs1*rs1) +A6/rs1-A7 + A8*(1.+log(rs1)));
        }
    }
}

void KS::set_Ecoul()
{
    
    for (int k = 0; k < num_points_mesh; k++){
        vector<double> I (num_points_mesh);
        for (int j = 0; j < num_points_mesh; j++)
            I[j] = (nOld[j]-n_ave)*abs(x_mesh_boundary[k] - x_mesh_boundary[j]);
        E[k] = -M_PI*e*e*(nOld[k]-n_ave)*tr.integra(I, x_mesh_boundary);
        //cout <<  tr.integra(I, x_mesh_boundary) << endl;
    }

    E_coul = tr.integra(E, x_mesh_boundary) ;
}

void KS::set_Eh()
{
    
    for (int k = 0; k < num_points_mesh; k++){
        vector<double> I (num_points_mesh);
        for (int j = 0; j < num_points_mesh; j++)
            I[j] = (nOld[j]-n_ave)*abs(x_mesh_boundary[k] - x_mesh_boundary[j]);
        E[k] = -2*M_PI*e*e*(nOld[k])*tr.integra(I, x_mesh_boundary);
        //cout <<  tr.integra(I, x_mesh_boundary) << endl;
    }

    E_h = tr.integra(E, x_mesh_boundary) ;
    
}

void KS::set_EH()
{
    for (int k = 0; k < num_points_mesh; k++){
        vector<double> I (num_points_mesh);
        for (int j = 0; j < num_points_mesh; j++)
            I[j] = (nOld[j])*abs(x_mesh_boundary[k] - x_mesh_boundary[j]);
        E[k] = -M_PI*e*e*(nOld[k])*tr.integra(I, x_mesh_boundary);
        //cout <<  tr.integra(I, x_mesh_boundary) << endl;
    }

    E_H = tr.integra(E,x_mesh_boundary);
}

void KS::set_Eb()
{
    for (int k = 0; k < num_points_mesh; k++){
        vector<double> I (num_points_mesh);
        for (int j = 0; j < num_points_mesh; j++)
            I[j] = abs(x_mesh_boundary[k] - x_mesh_boundary[j]);
        E[k] = 2*M_PI*e*e*n_ave*(nOld[k])*tr.integra(I, x_mesh_boundary);
        //cout <<  tr.integra(I, x_mesh_boundary) << endl;
    }
    E_b = tr.integra(E,x_mesh_boundary);
}



    
void KS::set_Exc()
{
    for (int k = 0; k <= dim; k++){
        if(k == 0){
            E[k] = 0;
        }else{
            if (k == (dim))
            {
                E[k] == 0;
            }else{
                E[k] =  Exc[k]*nOld[k];
            }
        }
         
    }
        
    E_xc = tr.integra(E,x_mesh_boundary);
}

void KS::set_EVxc()
{
    for (int k = 0; k <= dim; k++){
        if(k == 0){
            E[k] = 0;
        }else{
            if (k == (dim ))
            {
                E[k] == 0;
            }else{
                E[k] =  Vxc[k]*nOld[k];;
            }
        }
         
    }
    E_Vxc = tr.integra(E,x_mesh_boundary); 
}

void KS::set_Eext()
{
    for (int k = 0; k < num_points_mesh; k++)
        E[k] =  K*nOld[k]*(x_mesh_boundary[k]); 
    E_ext = tr.integra(E,x_mesh_boundary);
}



void KS::set_Ekin()
{
    E_kin = 0.;
    for (int k= 0; k<dim; k++)
    {   
        double I = 0;
        if (eigenvalues[k]<=eF ){
            for (int i = 0; i < dim; i++)
            {
               
               for (int j = 0; j <dim; j++){

                    I+=KS_matrix[ind(k,i)]*T[(ind(i,j))]*KS_matrix[ind(k,j)];
               }

            }
            //I = tr.integra(E,x_mesh_boundary);
            E_kin+=(eF-eigenvalues[k])*((eF-eigenvalues[k])*1/C-I);    
        } 
    }
    E_kin*=(1/(2*M_PI));
    
    
}


/*
void KS::set_Ekin()
{
    E_kin = 0.;
    for (int k= 0; k<dim; k++)
    {   
        double I = 0;
        if (eigenvalues[k]<=eF ){
            
                
            for (int i = 0; i < num_points_mesh; i++){
                if (i == 0){
                    E[i] = 0;
                }
                else{
                    if (i == num_points_mesh-1){
                        E[i] = 0;
                    }else{
                        E[i] = (KS_matrix[ind(k,i+1)] +KS_matrix[ind(k,i-1)]-2*KS_matrix[ind(k,i)]/(diff*diff)*KS_matrix[ind(k,i)]);
                    }
                }
            }
            I = tr.integra(E,x_mesh_boundary);
            E_kin+=(eF-eigenvalues[k])*((eF-eigenvalues[k])*1/C-I); 
        }
                                    
    } 
    E_kin*=(1/(2*M_PI));
}
*/

void KS::set_Et()
{
    E_t = 0.;
    for (int k= 0; k<dim; k++)
    {   
        if (eigenvalues[k]<=eF )
            E_t+=(eF*eF-eigenvalues[k]*eigenvalues[k]);     
    }
    
    E_t=E_t*(B/2.);
}

void KS::set_KS_potential()//settp potenziale a partire da densita' di input.
{
    //ofstream f {"DAT/pot.dat"};
   

    set_Vext();
    
    set_Vxc();
    
    for(r = 0;  r < dim; r++)
    {
        set_Vh();
        //Vnew[r] = tr.integra(Vh, x_mesh);
        Vnew[r] = tr.integra(Vh, x_mesh) + Vext[r] +Vxc[r];
        //Vnew[r] = tr.integra(Vh, x_mesh)+Vext[r];
        //Vnew[r] = Vext[r];
        //f << r << " " << tr.integra(Vh,x_mesh) << " " << Vext(r) << endl;
    };
    
}

void KS::set_KS_matrix()//setto la matrice di KS.
{
    for(int i=0; i<dim; ++i)
        for(int j=0; j<dim; ++j)
            KS_matrix[ind(i,j)] = 0;  
    for(int i=0; i<dim; ++i){
        KS_matrix[ind(i,i)]+= Vold[i];
        for(int j=0; j<dim; ++j)
            KS_matrix[ind(i,j)] += -(C/2.)*T[ind(i,j)];
            //KS_matrix[ind(i,j)]+= 0;
    }
    
}

void KS::norm()
{
    for(r = 0; r <dim; r++){//ciclo su tutti gli autovalori ottenuti dalla diagonalizzazione. Avremmo potuto ciclare solo su quelli con energia minore all'enrgia di fermi
        vector<double> Psi2 (dim);
        for (int k = 0; k < dim; k++)
            Psi2[k] = KS_matrix[ind(r,k)]*KS_matrix[ind(r,k)];

        double norm = tr.integra(Psi2, x_mesh);

        for (int k = 0; k < dim; k++) KS_matrix[ind(r,k)]/=sqrt(norm);   
    };
    
}


void KS::density()//setto la nuova densita'
{
    for(int j= + 0; j < dim; j++)
        nNew[j+1] = 0;
   
    for (r = 0;r <dim; r++){//ciclo su tutti gli autovalori ottenuti dalla diagonalizzazione. Avremmo potuto ciclare solo su quelli con energia minore all'enrgia di fermi
        double eFxy = B*(eF-eigenvalues[r]);
        if (eFxy > 0){//considero solo quelli con energia minore della energia di fermi
            for(int j= + 0; j < dim; j++){  //setto la densita' con gli autovettori ottenuti dalla diagonalizzazione relatvi a livelli occupati da elettroni in quanto sotto eF.
                nNew[j+1] += eFxy*pow(KS_matrix[ind(r,j)],2);
            } 
            nF = r+1; //conto il numero di livelli sotto l'energia di fermi: in concreto, le sottobande.
        } 
           
    };
    
    
}


void KS::mixing(vector<double>& vNew, vector<double>& vOld)//implemento mixing
{
    for (int i = 0; i < int(vOld.size()); i++)
        vOld[i] = (1-alpha)*vOld[i] + alpha*vNew[i];  

}

double KS::electron_counting(double ef_value)
{
    double sum =0;
    ofstream output_kFn {"DAT/kF.dat", ios::out};

    
    
    for (int n = 0; n < dim; n++)//ciclo su tutti gli autovalori ottenuti dalla diagonalizzazione
    {
        if (eigenvalues[n] <= ef_value){//consto solo quelli sotto l'energia di fermi data input
            double kF2 = (ef_value-eigenvalues[n]);
            
            //cout << "n = " << n << "; kFn = " << sqrt(kF2) << endl;
            output_kFn  << n << " " << kF2*B  << endl;
            sum+=kF2;
        }
        
    };
    sum*=B;
    
    
    return sum;
}

void KS::set_err()
{
    for(int k = 0; k < num_points_mesh; k++)
        Err[k] =  pow(nOld[k]-nNew[k],2);
}

void KS::eqKS()
{
        
    double errore;
    int conta = 0;
    set_init_density();
    
    do{
        
       
        

        set_KS_potential();
        
        //mixing(Vnew, Vold);
        for (int i = 0; i < int(Vold.size()); i++)
            Vold[i] = Vnew[i]; 
        
        //Vold = Vnew;
        
        set_KS_matrix();
        
        //print_file_KS_matrix(0);
        diagonalize(dim, KS_matrix, eigenvalues, lwork, work, info);
        //print_file_KS_matrix(1);
        
        norm();
        
        

        double ef_inf = eigenvalues[0];
        double ef_sup = eigenvalues[dim-1];

        

        eF = bi.esegui(ef_inf, ef_sup, toll,fsigma);

        print_file_eigenvalues();

        

        
       

        m_sigma = fsigma(eF)+sigma;
        m_n_ave = m_sigma/L;
        
        density();

        set_err();
        errore = tr.integra(Err,x_mesh_boundary);

        mixing(nNew, nOld);
        
        //nOld = nNew;
         
        print_file_density();
       
        set_exc();
        
        set_Ekin();//CALCOLO ENERGIA CINETICA: CONTO ESPLICITO
        set_Ecoul();//CALCOLO TERMINE COULOMBIANO
        set_Eext();//ENERGIA ESTERNA
        set_Exc();//ENERGIA SCAMBIO E CORRELAZIONE

        set_Et();//SOMMA AUTOVALORI
        set_EVxc();//ENERGIA TCAMBIO E CORRELAZIONE TERMINE DERIVATA
        set_Eh();//TERMINE HYBRID: COULOMBIANO + INTEGRAZIONE RISPETTO DENSITA'
        set_EH();//TERMINE DI HARTREE IN SENSO STRETTO
        set_Eb();//TERMINE DI BACKGROUND DI INTERAZIONE CON GLI ELETTRONI
        

        Enew = E_kin + E_coul + E_ext + E_xc;
        //Enew =  +E_t  -E_H;
        //Enew =  -E_h;
        //E_kin = E_t-E_h-E_H -E_ext-E_Vxc;

        
        if (conta %20==0){
            
            cout << "$$$$$$$$$$$$$$$" << endl;
            cout << "ef_inf = " << ef_inf<< endl;
            cout << "ef_sup= "<<ef_sup << endl;
            cout << "fsigma(ef_inf) = " << fsigma(ef_inf) << endl;
            cout << "fsigma(ef_sup) = " << fsigma(ef_sup) << endl;
            cout << "ef = " << eF << endl;
            cout << "sigma = " << fsigma(eF)+sigma << endl;
            cout << "nF = " << nF << endl; 
            cout << "AUTOVALORI: ";
            for (int i = 0; i < nF; i++) cout << eigenvalues[i] << " ";
            cout << endl;
            cout << "n_ave = " << tr.integra(nOld, x_mesh_boundary)/L << endl;
            cout << "sigma = " << tr.integra(nOld, x_mesh_boundary) << endl;
            cout << "Ekin = " << E_kin << endl;
            cout << "Ecoul =" << E_coul << endl;
            cout << "Eext = " << E_ext << endl;
            cout << "Exc = " << E_xc << endl;
            cout << "ENERGIA 1 = " << E_kin + E_coul + E_ext + E_xc << endl;
            cout << "Et = " << E_t << endl;
            cout << "Eh =" << E_h << endl;  
            cout << "EVxc = " << E_Vxc << endl;
            cout << "Eb = " << E_b << endl;
            cout << "EH = " << E_H << endl; 
            cout << "Ebb = " << E_bb << endl;
            cout << "err = " << errore << endl;
        
        
            cout << "$$$$$$$$$$$$$$$" << endl;
        
            cout << "$$$$$$$$$$$$$$$" << endl;
        }
        
        conta++;
    }while(errore > prec  );
    

   
    
}




//METODI PRINT

void KS::print_file_KS_matrix(int a)
{
    ofstream output_KS_matrix {"DAT/KS_matrix"+to_string(a)+".dat"};
    
    for(int i=0; i<dim; ++i){
        for(int j=0; j<dim; ++j)
            output_KS_matrix << KS_matrix[ind(i,j)] << " ";
            //KS_matrix[ind(i,j)]+= 0;
        output_KS_matrix << endl;
    }
}
void KS::print_file_sigma_alpha()
{
    ostringstream oss;
    oss << setprecision(2)<< K; 
    string a_K = oss.str();
    oss.str("");
    string title = "DAT/sigma_alpha_K="+a_K+".dat";
    
    fstream output_sigma_alpha {title.c_str(), ios::app};
    for (int n = 0; n < dim; n++)//ciclo su tutti gli autovalori ottenuti dalla diagonalizzazione
    {
        if (eigenvalues[n] <= eF){//consto solo quelli sotto l'energia di fermi data input
            double kF2 = (eF-eigenvalues[n]);
            
            //cout << "n = " << n << "; kFn = " << sqrt(kF2) << endl;
            if (n == 0)
                sigma_alpha1 = kF2*B;
            if (n == 1)
                sigma_alpha2 = kF2*B;
            if(n == 2)
                sigma_alpha3 = kF2*B ;
        }
        
    };

    output_sigma_alpha << sigma << " " << nF << " " << sigma_alpha1 << " " << sigma_alpha2 << " " << sigma_alpha3 << endl;
}


void KS::print_file_T()//STAMPA SU FILE MATRICE CINETICA
{
    ofstream T_out {"DAT/T.dat"};
    for(int i = 0; i <dim; i++){
        for(int j = 0; j < dim; j++)
            T_out << T[ind(i,j)] << " ";
        T_out << endl;
    }
    T_out.close();
}

void KS::print_file_eigenvalues_result()//STAMPA SU FILE ATUOVALORI
{
    ofstream output_eigenvalues {"DAT/eigenvalues_result.dat", ios::app};
    output_eigenvalues << n_ave << " " << eF << " ";
    for (int k = 0; k < dim; k++){
        output_eigenvalues <<  eigenvalues[k] << " ";
    }

    output_eigenvalues << endl;
    output_eigenvalues.close();
        
}

void KS::print_file_eigenvalues()//STAMPA SU FILE ATUOVALORI
{
    ofstream output_eigenvalues {"DAT/eigenvalues.dat"};
    output_eigenvalues << "eF = " << eF << endl;
    for (int i = 0; i < dim; i++){
            output_eigenvalues <<  eigenvalues[i] << endl;
            
    }
    
    output_eigenvalues.close();
        
}

void KS::print_file_density()//STAMPA SU FILE DENSITA'
{
    ostringstream oss;
    oss << setprecision(2)<< K; 
    string a_K = oss.str();
    oss.str("");
    oss << setprecision(2) << sigma;
    string a_sigma = oss.str();
    oss.str("");
    string title = "DAT/density_K="+a_K+"_sigma="+a_sigma+".dat";
    
    fstream output_density {title.c_str(), ios::out};
    //ofstream output_density {"DAT/density.dat"};
    for(int k = 0; k < dim+2; k++){
        
        output_density << x_mesh_boundary[k] << " " <<  nOld[k] << endl;
        
    }
    
    output_density.close();
    
}

void KS::print_file_energy()//STAMPA SU FILE ENERGIA
{
    //string title = "DAT/energy_K="+to_string(int(K))+"N="+to_string(num_part)+".dat";
    
    ofstream output_energy {"DAT/energy.dat", ios::app};
    output_energy << sigma << " " << Enew << " " << E_kin << " " << E_coul << " " << E_ext << " " << E_xc<< " " << E_Vxc << " " << E_H << " " << E_t << " " << E_b << " " << E_bb << endl;
    output_energy.close();
}

void KS::print_file_subbands()
{
    ofstream output_subbands {"DAT/subbands.dat", ios::app};
    output_subbands << sigma << " " << nF << endl;
    output_subbands.close();
}

void KS::print_file_eigenfunctions(int a)
{
    ostringstream oss;
    oss << setprecision(2)<< K; 
    string a_K = oss.str();
    oss.str("");
    oss << setprecision(2) << sigma;
    string a_sigma = oss.str();
    oss.str("");
    for (int j = 1; j <= a; j++)
    {
        string title = "DAT/eigenfunction_K="+a_K+"_sigma="+a_sigma+"_"+to_string(j)+".dat";
        fstream f {title.c_str(), ios::out};
        for (int i = 0; i < dim; i++)
            f << x_mesh[i] << " " << pow(KS_matrix[ind(j-1,i)],2) << endl;
    }
}