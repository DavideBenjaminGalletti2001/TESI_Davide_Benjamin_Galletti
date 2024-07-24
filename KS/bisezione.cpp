#include "./LIB/bisezione.h"

//COSTRUTTORE
Bisezione::Bisezione(int iter_max , double toll )
{

    m_iter_max = iter_max;
    m_toll = toll;
}



//FUNZIONE SEGNO
int Bisezione::sng(double &x)
{
    int ris = x > 0. ? 1 :( x < 0 ? -1 : 0);
    return ris;
}

//CERCA ZERI
double Bisezione::esegui(std::function <double (double)> f)
{
    m_iter = 0;
    double m_e = std::abs(m_b - m_a);
    
    double fa = f(m_a);
    double fb = f(m_b);
    
    double fx ;
    double a = m_a;
    double b = m_b;

    while (m_e > m_toll)
    {
        if (m_iter > m_iter_max ) throw std::runtime_error("Error. Number of iterations overcames.");
        m_x = (a + b)/2.;
        fx = f(m_x);
        if(sng(fa)*sng(fx) < 0)
        {
            b = m_x;
            fb = fx;
        }
        else
        {
            if (sng(fb)*sng(fx)<=0)
            {
                a = m_x;
                fa = fx;
            }
            else
            {
                return m_x;
            }
        }
        m_iter++;
        m_e = std::abs(b - a);
    }
    m_e=m_e/2.;
    return m_x;   
}


//CERCA ZERI
double Bisezione::esegui(double & a, double & b, double & toll, std::function <double (double)> f)
{
    m_iter = 0;

    
    m_a = a;
    m_b = b;
    m_e = std::abs(m_b - m_a);
    double fa = f(m_a);
    double fb = f(m_b);
    

    if (sng(fa)*sng(fb)>0)  throw std::runtime_error("Error. The extremes don't satisfay the hypotesis!");
    
    double fx ;
    m_toll = toll;
    m_iter_max = 100000;
    while (m_e > m_toll)
    {
        if (m_iter > m_iter_max ) throw std::runtime_error("Error. Number of iterations overcames.");
        m_x = (a + b)/2.;
        fx = f(m_x);
        
        if(sng(fa)*sng(fx) < 0)
        {
            b = m_x;
            fb = fx;
            
        }
        else
        {
            if (sng(fb)*sng(fx)<=0)
            {
                a = m_x;
                fa = fx;
                
            }
            else
            {
                
                return m_x;
            }
        }
        m_iter++;
        m_e = std::abs(b - a);
        
    }
    m_e=m_e/2.;
    return m_x;   
}
void Bisezione::print(void)
{
    std::cout << "ZERO = " << m_x << " ERRORE = " << m_e << std::endl;
}