//LSN exercise 02.1
//Computing a 1D integral via Monte Carlo:
//𝐼=∫𝜋/2cos(𝜋𝑥/2)𝑑𝑥=1


#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"

int main ()
{
    const unsigned int M = 1E5 ; //Number of throws
    const unsigned int N = 1E2 ;  //Number of blocks
    unsigned int L = M / N ; //Number of throws in each block
    
    double x = 0. ;  //Random numbers generated
    double f = 0. ;  //Evaluation of the function
    
    //Average of every block, Square of the average in every block, Standard Deviation
    double A = 0., A2 = 0., sigma = 0. ;
    
    Random gen ("Primes", "seed.in") ;  //Random numbers generator

    
    //1. Sampling a uniform distribution in  [0,1]
    
    //Creating output file
    std::ofstream WriteUnif ;
    WriteUnif.open ( "Integral_Unif.dat" ) ;
    
    if ( !WriteUnif.is_open() )
    {
        std::cerr << "attention! Can't create output file for uniform distribution!" ;
        return 1 ;
    }
    
     //Generating the random numbers and evaluating the function, calculating the average with the blocks method and writings the results on the file
    for ( int j = 0; j < N; j ++ )
    {
        for ( int i = 0; i < L; i ++ )
        {
            x = gen.Rannyu() ;
            f += (M_PI/2.) * cos(M_PI * x / 2.) ;
        }
        f /= L ;
        A += f ;
        A2 += f * f ;
        sigma = ( A2 / (double)(j+1.) ) - pow ( A / (double)(j+1.), 2) ;
        if (j==0)
            sigma = 0;
        else
            sigma = std::sqrt (sigma / (double)j ) ;
        
        WriteUnif << j << " " << A/(j+1.) << " " << sigma << std::endl ;
    
        f = 0. ;
    }
    
    WriteUnif.close() ;
    
    
    //2. Using importance sampling
    
    
    //Creating output file
    std::ofstream WriteImp ;
    WriteImp.open ( "Integral_Imp.dat" ) ;
       
       if ( !WriteImp.is_open() )
       {
           std::cerr << "attention! Can't create output file for importance sampling!" ;
           return 1 ;
       }
    
    //Re-initialising
    A = 0;
    A2 = 0;

     //Generating the random numbers following the new dostribution and evaluating the new function, calculating the average with the blocks method and writing the results on the file
    for ( int j = 0; j < N; j ++ )
    {
    
        for ( int i = 0; i < L; i ++ )
        {
            x = 1 - sqrt(gen.Rannyu() ) ;
            f += (M_PI/2.) * cos(M_PI * x / 2.) / ( 2. * (1. - x) ) ;
        }
        f /= L ;
        A += f ;
        A2 += f * f ;
        sigma = ( A2 / (double)(j+1.) ) - pow ( A / (double)(j+1.), 2) ;
        
        if (j==0)
            sigma = 0;
        else
            sigma = std::sqrt (sigma / (double)j ) ;
  
            WriteImp << j << " " << A/(j+1.) << " " << sigma << std::endl ;
       
        f = 0. ;
    }
    
    WriteImp.close() ;
    
    return 0;
}
