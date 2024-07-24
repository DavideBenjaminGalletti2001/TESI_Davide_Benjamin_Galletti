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
#include "solutore.h" //CLASSE MADRE ASTRATTA

#ifndef __bisezione__
#define __bisezione__

class Bisezione{
    public:
        Bisezione(double a = 0, double b = 1, unsigned int iter = 10000, double toll = 0.001);
        double esegui(std::function <double (double)> f);
        double esegui(double  a, double  b, double & toll, std::function <double (double)> f);
        void check_interval(double & a, double & b);
        int sng (double & x);
        void print(void);
    private:
        double m_a, m_b;
        double m_e ;
        unsigned int m_iter;
        unsigned int m_iter_max;
        double m_toll;
        double m_x;
};
#endif