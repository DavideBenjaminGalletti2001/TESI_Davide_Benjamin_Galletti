#ifndef __trapezi__
#define __trapezi__
#include "integrali.h"
class trapezi: public integral{
    public:
        trapezi(): integral(){;};
        double integra(std::function <double (int)> f, vector<double> &x);
        double integra_runtime(double prec, std::function <double (int)> f, vector<double> &x);
        double integra(vector<double>& f, vector<double> &x);
        double integra_runtime(double prec, vector<double>& f, vector<double> &x);
    private:
};
#endif