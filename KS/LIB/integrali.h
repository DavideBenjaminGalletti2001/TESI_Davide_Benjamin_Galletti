#ifndef __integrali_h__
#define __integrali_h__

#include <cmath>
#include <vector>
#include<iostream>
#include <iomanip>
#include <cstdlib>
#include <functional>

using namespace std;

class integral{
    public:
        integral();
        
        virtual double integra(std::function<double (int)> f, vector<double> &x) = 0;
        virtual double integra_runtime(double prec, std::function<double (int)> f, vector<double> &x)  = 0;
        virtual double integra(vector<double>& f, vector<double> &x) = 0;
        virtual double integra_runtime(double prec, vector<double>& f, vector<double> &x)  = 0;
        
    

    protected:
        void checkInterval(double a, double b);
        int m_nstep;
        double m_a, m_b;
        double m_sum, m_integral, m_h;
        int m_sign;
        double m_prec;
};




#endif