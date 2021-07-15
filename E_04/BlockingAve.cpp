
#include "BlockingAve.hpp"

Blocking::Blocking ()
{
    m_counter = 0 ;
    m_val = 0. ;
    m_A = 0. ;
    m_A2 = 0. ;
    m_sigma = 0. ;
    m_ave = 0. ;
    m_N = 0 ;
    m_L = 0 ;
}

Blocking::Blocking (int M, int N)
{
    m_counter = 0 ;
    m_val = 0. ;
    m_A = 0. ;
    m_A2 = 0. ;
    m_sigma = 0. ;
    m_ave = 0. ;
    m_N = N ;
    m_L = M/N ;
    if (M % N != 0)
    std::cerr << "WARNING: M % N != 0. Please change the values of M and N so that M is divisible by N." << std::endl << std::endl ;
    }


void Blocking::Add (double n, std::ofstream & Write)
{
    m_val += n ;
    m_counter ++ ;
    
    int j = m_counter / m_L ;
    
    if (m_counter % m_L == 0)
    {
        m_val /= (double)m_L ;
        m_A += m_val ;
        m_A2 += m_val * m_val ;
        m_sigma = ( m_A2 / (double)(j) ) - pow ( m_A / (double)(j), 2) ;
        
        if (j==1)
            m_sigma = 0;
        else
            m_sigma = sqrt (m_sigma / (double)(j-1) ) ;
       
        Write << j-1 << " " << std::setprecision(9) << m_A / (double)(j) << " " << m_sigma << std::endl ;
        
        m_ave = m_A / (double)(j) ;
        m_val = 0. ;
    }
}

void Blocking::Writefinal (double n, double nblk, std::ofstream  & Write)
{
    m_val += n ;
    m_counter ++ ;
    
    int j = m_counter / m_L ;
    
    if (m_counter % m_L == 0)
    {
        m_val /= (double)m_L ;
        m_A += m_val ;
        m_A2 += m_val * m_val ;
        m_sigma = ( m_A2 / (double)(j) ) - pow ( m_A / (double)(j), 2) ;
        
        if (j==1)
            m_sigma = 0;
        else
            m_sigma = sqrt (m_sigma / (double)(j-1) ) ;
        
        if (j == nblk)
            Write << m_A / (double)(j) << " " << m_sigma << std::endl ;
        
        m_ave = m_A / (double)(j) ;
        m_val = 0. ;
    }
}





void Blocking::Reset ()
{
    m_counter = 0 ;
    m_val = 0. ;
    m_A = 0. ;
    m_A2 = 0. ;
    m_sigma = 0. ;
    m_ave = 0. ;
}


void Blocking::Set_parameters (int M, int N)
{
    m_N = N ;
    m_L = M/N ;
    if (M % N != 0)
        std::cerr << "WARNING: M % N != 0. Please change the values of M and N so that M is divisible by N." << std::endl << std::endl ;
}
