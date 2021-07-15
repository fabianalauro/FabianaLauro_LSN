//
//  PlainVanilla.cpp
//  
//
//  Created by Fabiana on 06/04/21.
//

#include "PlainVanilla.hpp"



PlainVanilla::PlainVanilla ()
{
    m_gen = NULL ;
    m_t0 = 0. ;
    m_t1 = 0. ;
    m_S0 = 0. ;
    m_S1 = 0. ;
    m_K = 0. ;
    m_r = 0. ;
    m_sigma = 0. ;
}

PlainVanilla::PlainVanilla (double t0, double t1, double S0, double K, double r, double sigma, Random * gen)
{
    m_gen = gen ;
    m_t0 = t0 ;
    m_t1 = t1 ;
    m_S0 = S0 ;
    m_S1 = 0. ;
    m_K = K ;
    m_r = r ;
    m_sigma = sigma ;
}
   
void PlainVanilla::setRandom (Random * gen)
{
    m_gen = gen ;
}

void PlainVanilla::setParameters (double t0, double t1, double S0, double K, double r, double sigma)
{
    m_t0 = t0 ;
    m_t1 = t1 ;
    m_S0 = S0 ;
    m_K = K ;
    m_r = r ;
    m_sigma = sigma ;
}

void PlainVanilla::setS1_dir ()
{
    m_S1 = m_S0 * exp ( (m_r - m_sigma*m_sigma/2.) * (m_t1 - m_t0) + m_sigma * m_gen->Gauss(0., m_t1 - m_t0) ) ;
}

void PlainVanilla::setS1_discr (unsigned int n_steps)
{
    double S = m_S0 ;
    double t_incr = (m_t1 - m_t0) / n_steps ;
    for (int i = 0; i < n_steps; i++)
    {
        S *= exp ( (m_r - m_sigma*m_sigma/2.) * (t_incr) + m_sigma * m_gen->Gauss(0., 1.) * sqrt(t_incr) ) ;
    }
    m_S1 = S ;
}

double PlainVanilla::S1 ()
{
    return m_S1 ;
}
   
double PlainVanilla::Call_Profit ()
{
    return exp (- m_r * m_t1) * fmax ( 0., m_S1 - m_K ) ;
}

double PlainVanilla::Put_Profit ()
{
    return exp (- m_r * m_t1) * fmax ( 0., m_K - m_S1 ) ;
}
