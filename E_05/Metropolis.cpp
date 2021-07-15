//

#include "Metropolis.hpp"

Metropolis3D::Metropolis3D ()
{
    m_x = 0. ;
    m_y = 0. ;
    m_z = 0. ;
    m_alpha = 0. ;
    m_gen = NULL ;
    m_prob = NULL ;
}

Metropolis3D::Metropolis3D (double x, double y, double z)
{
    m_x = x ;
    m_y = y ;
    m_z = z ;
    m_alpha = 0. ;
    m_gen = NULL ;
    m_prob = NULL ;
}


Metropolis3D::Metropolis3D (double x, double y, double z, Random* gen, Basic_PDF * prob)
{
    m_x = x ;
    m_y = y ;
    m_z = z ;
    m_alpha = 0. ;
    m_gen = gen ;
    m_prob = prob ;
}

double Metropolis3D::r () const
{
    return sqrt (m_x*m_x + m_y*m_y + m_z*m_z) ;
}

double Metropolis3D::alpha () const
{
    return m_alpha ;
}

 void Metropolis3D::Set_x (double n)
{
    m_x = n ;
    return ;
}

void Metropolis3D::Set_y (double n)
{
    m_y = n ;
    return ;
}

void Metropolis3D::Set_z (double n)
{
    m_z = n ;
    return ;
}

void Metropolis3D::Set_coord (double x, double y, double z)
{
    m_x = x ;
    m_y = y ;
    m_z = z ;
    return ;
}

void Metropolis3D::Set_PDF (Basic_PDF * prob) 
{
    m_prob = prob ;
    return ;
}





void Metropolis3D::MRT2_Unif (double max)
{
    double xnew = m_x + m_gen->Rannyu(-max, max) ;
    double ynew = m_y + m_gen->Rannyu(-max, max) ;
    double znew = m_z + m_gen->Rannyu(-max, max) ;
    m_alpha = fmin ( 1, m_prob->Probability(xnew, ynew, znew) / m_prob->Probability(m_x, m_y, m_z) ) ;
    double r = m_gen->Rannyu() ;
    if ( r < m_alpha || r ==m_alpha )
    {
        m_x = xnew ;
        m_y = ynew ;
        m_z = znew ;
    }
    return ;
}

void Metropolis3D::MRT2_Gauss (double sigma)
{
    double xnew = m_x + m_gen->Gauss(0, sigma) ;
    double ynew = m_y + m_gen->Gauss(0, sigma) ;
    double znew = m_z + m_gen->Gauss(0, sigma) ;
    m_alpha = fmin ( 1, m_prob->Probability(xnew, ynew, znew) / m_prob->Probability(m_x, m_y, m_z) ) ;
    double r = m_gen->Rannyu() ;
    if ( r < m_alpha || r == m_alpha )
    {
        m_x = xnew ;
        m_y = ynew ;
        m_z = znew ;
    }
    return ;
}


