//Classes to implement the metropolis algorithm, in various dimensions

#ifndef Metropolis_hpp
#define Metropolis_hpp

#include <cmath>
#include "random.h"
#include "Probability.hpp"  //Where the PDF is written (a class with that specific PDF needs to be created, as a derived class of Basic_PDF)



class Metropolis1D
{
private:
    
    double m_x ;
    double m_alpha ;
    Random * m_gen ;
    Basic_PDF * m_prob ;
    
protected:
    
    
public:
    
    Metropolis1D () ;
    Metropolis1D (double x) ;
    Metropolis1D (double x, Random * gen, Basic_PDF * prob) ;
    
    double x () const {return m_x;} ;
    double alpha () const ;
    
    void Set_x (double n) ;
    void Set_PDF (Basic_PDF * prob) ;
    
    void MRT2_Unif (double max) ; //Implementation of Metropolis algorithm using uniform distribution
    void MRT2_Gauss (double sigma) ; //Implementation of Metropolis algorithm using gaussian distribution
    
    
} ;



class Metropolis3D
{
private:
    
    double m_x, m_y, m_z ;
    double m_alpha ;
    Random * m_gen ;
    Basic_PDF * m_prob ;
    
protected:
    
    
public:
    
    Metropolis3D () ;
    Metropolis3D (double x, double y, double z) ;
    Metropolis3D (double x, double y, double z, Random * gen, Basic_PDF * prob) ;
    
    double x () const {return m_x;} ;
    double y () const {return m_y;} ;
    double z () const {return m_z;} ;
    double r () const ;
    double alpha () const ;
    
    void Set_x (double n) ;
    void Set_y (double n) ;
    void Set_z (double n) ;
    void Set_coord (double x, double y, double z) ;
    void Set_PDF (Basic_PDF * prob) ;
    
    void MRT2_Unif (double max) ; //Implementation of Metropolis algorithm using uniform distribution
    void MRT2_Gauss (double sigma) ; //Implementation of Metropolis algorithm using normal distribution
    
    
} ;
















#endif /* Metropolis_hpp */
