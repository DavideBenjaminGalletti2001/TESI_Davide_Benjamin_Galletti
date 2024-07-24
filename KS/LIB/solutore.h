#ifndef __solutore__
#define __solutore__
#include <functional>
#include <iostream>
class Solutore{
    public:
        Solutore();
        virtual double esegui(std::function <double (double)> f) = 0;
        virtual double esegui(double & a, double & b, double & toll, std::function <double (double)> f) = 0;
    protected:
};
#endif
