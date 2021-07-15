//Class to compute Blocking Average (without using std::vectors) and writing the results on a file


#ifndef BlockingAve_hpp
#define BlockingAve_hpp

#include <iostream>
#include <fstream>
#include <cmath> //pow, sqrt

class Blocking
{
private:
    
    int m_counter ;
    double m_val, m_A, m_A2, m_sigma, m_ave ;
    int m_N, m_L ;
  
protected:
    
public:
    
    Blocking () ;
    Blocking (int M, int N) ;

    
    void Set_parameters (int M, int N) ;
    double ave () const {return m_ave;};
    double err () const {return m_sigma;};
    
    void Reset () ;
    void AddProg (double n, std::ofstream  & Write) ; //Performs blocking average and prints on a file the progressive averages
    void AddFinal (double n, double nblk, std::ofstream  & Write) ; //Performs blocking average and prints on a file only the final value of the average
    void Add (double n) ; //Performs blocking average 
};





#endif /* BlockingAve_hpp */
