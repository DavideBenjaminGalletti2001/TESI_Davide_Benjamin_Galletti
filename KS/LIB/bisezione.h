#include "solutore.h"
#ifndef __bisezione__
#define __bisezione__
#include <iostream>
#include <fstream>
class Bisezione{
    public:
        Bisezione(int iter = 10000, double toll = 0.001);
        double esegui(std::function <double (double)> f);
        double esegui(double & a, double & b, double & toll, std::function <double (double)> f);
        
        int sng (double & x);
        void print(void);
    private:
        double m_a, m_b;
        double m_e ;
        int m_iter;
        int m_iter_max;
        double m_toll;
        double m_x;
};
#endif