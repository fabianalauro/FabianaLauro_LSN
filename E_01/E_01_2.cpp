//LSN Exercise 02.1
//Checking the Central Limit Theorem


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"


int main ()
{
    Random gen ("Primes", "seed.in") ;
    
    const unsigned int M = 1E4 ;  //Total number of realisations of S_N
    double x = 0. , y = 0. , z = 0. ; //The random numbers generated for each distribution
    
    double N[] = {1., 2., 10., 100.} ;
    
    //Creating the output files for each distributioin
    std::ofstream WriteU ;
    WriteU.open ("Unif.dat") ;
    if (!WriteU.is_open())
    {
        std::cerr << "Warning! Can't open Unif output file" ;
        return 1 ;
    }
    
    std::ofstream WriteE ;
    WriteE.open ("Exp.dat") ;
    if (!WriteE.is_open())
    {
        std::cerr << "Warning! Can't open Exp output file" ;
        return 1 ;
    }
       
    std::ofstream WriteL;
    WriteL.open ("Lor.dat") ;
    if (!WriteL.is_open())
    {
        std::cerr << "Warning! Can't open Lor output file" ;
        return 1 ;
    }
    
    //Generating the random numbers and writing them on the files
    for (int n = 0; n < 4; n++)
    {
        for (int i = 0; i < M; i++)
        {
            for (int j = 0; j < N[n]; j++)
            {
                x += gen.Rannyu (1., 6.) ;
                y += gen.Exp (1.) ;
                z += gen.Lorentz (0., 1.) ;
            }
            WriteU << x / N[n] << std::endl ;
            WriteE << y / N[n] << std::endl ;
            WriteL << z / N[n] << std::endl ;
            //Re-initialising
            x = 0. ;
            y = 0. ;
            z = 0. ;
        }
    }
    
    WriteU.close() ;
    WriteE.close() ;
    WriteL.close() ;
    
    return 0 ;
}
