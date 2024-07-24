#include "./LIB/integrali.h"

integral::integral(){
    m_sum = 0;
    m_integral = 0;
};

void integral::checkInterval(double a, double b){
    m_a = std::min(a,b);
    m_b = std::max(a,b);
    if (a>b)m_sign = -1;
    else m_sign = 1;
}
