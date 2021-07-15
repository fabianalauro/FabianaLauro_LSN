//Class to compute Blocking Average (without using std::vectors) and writing the results on a file


#ifndef BlockingAve_hpp
#define BlockingAve_hpp

#include <iostream>
#include <fstream> //ofstream
#include <cmath> //pow, sqrt
#include <iomanip> // setprecision

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
    
    void Reset () ;
    void Add (double n, std::ofstream  & Write) ; //Performs blocking average and prints on a file the progressive averages
    void Writefinal (double n, double nblk, std::ofstream  & Write) ; //Performs blocking average and prints on a file only the final value of the average
    
};


#endif /* BlockingAve_hpp */
