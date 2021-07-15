
#include "Probability.hpp"


double Hydr100::Probability (double x, double y, double z) 
{
    return 1./M_PI * exp ( -2. * sqrt (x*x + y*y + z*z) ) ;
}

double Hydr210::Probability (double x, double y, double z)
{
    return 1./(32 * M_PI) * z*z * exp ( - sqrt (x*x + y*y + z*z) ) ;
}

Variational::Variational ()
{
    m_sigma = 0. ;
    m_mu = 0. ;
}

Variational::Variational (double sigma, double mu)
{
    m_sigma = sigma ;
    m_mu = mu ;
}

void Variational::set_parameters (double sigma, double mu)
{
    m_sigma = sigma ;
    m_mu = mu ;
}

double Variational::Probability (double x, double y, double z)
{
    return pow ((exp(-pow(x-m_mu, 2) / (2*m_sigma*m_sigma) ) + exp (- pow(x+m_mu, 2) /(2* m_sigma*m_sigma ))),2) ;
}
