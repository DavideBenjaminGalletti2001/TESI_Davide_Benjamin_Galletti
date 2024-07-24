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
#include "./LIB/trapezi.h"

double trapezi::integra(std::function <double (int)> f, vector<double> &x)
{
    vector <double> m_delta_x;
    for (int i = 0; i < int(x.size())-1; i++) m_delta_x.push_back(x[i+1] - x[i]);
    m_sum = 0.;
    for (int i = 0; i < int(x.size())-1; i++)
    {
        
        m_sum += (f(i+1)+ f(i))*m_delta_x[i];
                 
    }
    

    m_integral = m_sum/2.;
    return m_integral;
}

double trapezi::integra(vector<double> & f, vector<double> &x)
{
    vector <double> m_delta_x;
    for (int i = 0; i < int(x.size())-1; i++) m_delta_x.push_back(x[i+1] - x[i]);
    m_sum = 0.;
    for (int i = 0; i < int(x.size())-1; i++)
    {
        
        m_sum += (f[i+1]+ f[i])*m_delta_x[i];
                 
    }
    

    m_integral = m_sum/2.;
    return m_integral;
}

double trapezi::integra_runtime(double prec, std::function <double (int)> f, vector<double> &x)
{
    if (prec <= 0) {std::cerr << "errore: il numero degli step è negativo" << std::endl;exit(-1);};

    m_prec = prec;

    double I1, I2;

    double e;

    m_nstep = 2;

    do{
        I1 = integra(f, x);
        I2 = integra(f, x);

        e = 4./3.*fabs(I1-I2);

    }while(e > m_prec);

    m_integral = I1;
    return e;
}

double trapezi::integra_runtime(double prec, vector<double> & f, vector<double> &x)
{
    if (prec <= 0) {std::cerr << "errore: il numero degli step è negativo" << std::endl;exit(-1);};

    m_prec = prec;

    double I1, I2;

    double e;

    m_nstep = 2;

    do{
        I1 = integra(f, x);
        I2 = integra(f, x);

        e = 4./3.*fabs(I1-I2);

    }while(e > m_prec);

    m_integral = I1;
    return e;
}