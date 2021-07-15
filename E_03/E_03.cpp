//LSN Exercise 03

//Make four pictures for the estimation of the European call-option prices,  𝐶[𝑆(0),0]  (direct and discretized), and put-option prices,  𝑃[𝑆(0),0]  (direct and discretized), with their uncertainties with a large number of asset prices at time  𝑡=𝑇 , say  𝑀≥10^4 , as a function of the number of blocks,  𝑁 . As usual, in computing statistical uncertainties, use data blocking.


#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "PlainVanilla.hpp"



int main ()
{
    //Constants for the blocking method
    const unsigned int M = 1E5 ;  //Number of throws
    const unsigned int N = 100 ;   //Number of blocks
    unsigned int L = M / N ;   //Number of throws in each block
    
    //Variables to average inside the block
    double Call = 0., Put = 0. ;
    
    //Average of every block, Square of the average in every block, Standard Deviation
    double A_Call = 0., A2_Call = 0., sigma_Call = 0. ;
    double A_Put = 0., A2_Put = 0., sigma_Put = 0. ;
    
    //Parameters
    const double t0 = 0. ;
    const double S0 = 100. ;
    const double T = 1 ;
    const double K = 100 ;
    const double r = 0.1 ;
    const double sigma = 0.25 ;
    
    //Random numbers generator
    Random gen ("primes32001.in", "seed.in") ;
    Random * p_gen = &gen ;
    
    
    //Class PlainVanilla
    PlainVanilla option (t0, T, S0, K, r, sigma, p_gen) ;
    
    
    //1. Direct method
    
    //Creating output file
    std::ofstream WriteDir ;
    WriteDir.open ( "Results_dir.dat" ) ;
    
    if ( !WriteDir.is_open() )
    {
        std::cerr << "attention! Can't create \"direct\" output file!" ;
        return 1 ;
    }
    
    
    //Simulation using direct method (PUT + CALL)
    for (int j = 0; j < N; j ++)
    {
        for ( int i = 0; i < L; i ++ )
        {
            option.setS1_dir () ;
            Call += option.Call_Profit () ;
            Put += option.Put_Profit () ;
        }
        Call /= L ;
        Put /= L ;
        
        //Averaging with blocking method
        A_Call += Call ;
        A2_Call += Call*Call ;
        sigma_Call = ( A2_Call / (double)(j+1.) ) - pow ( A_Call / (double)(j+1.), 2) ;
        if (j==0)
            sigma_Call = 0;
        else
            sigma_Call = std::sqrt (sigma_Call / (double)j ) ;
        
        A_Put += Put ;
        A2_Put += Put*Put ;
        sigma_Put = ( A2_Put / (double)(j+1.) ) - pow ( A_Put / (double)(j+1.), 2) ;
        if (j==0)
            sigma_Put = 0;
        else
            sigma_Put = std::sqrt (sigma_Put / (double)j ) ;
          
        //Writing on the file
        WriteDir << j << " " << A_Call/(j+1.) << " " << sigma_Call << " " << A_Put/(j+1.) << " " << sigma_Put << std::endl ;
          
    }
    
    WriteDir.close () ;
    
    
    //Re-initialising
    Call = 0.;
    Put = 0. ;
    A_Call = 0.;
    A2_Call = 0.;
    sigma_Call = 0. ;
    A_Put = 0.;
    A2_Put = 0.;
    sigma_Put = 0. ;
    
    //2. Discretised method

    //Creating output file
    std::ofstream WriteDiscr ;
    WriteDiscr.open ( "Results_discr.dat" ) ;
      
    if ( !WriteDiscr.is_open() )
    {
        std::cerr << "attention! Can't create \"discretised\" output file!" ;
        return 1 ;
    }
      
      
    //Simulation using discretsed method (PUT + CALL)
    for (int j = 0; j < N; j ++)
    {
        for ( int i = 0; i < L; i ++ )
        {
            option.setS1_discr (100) ;
            Call += option.Call_Profit () ;
            Put += option.Put_Profit () ;
        }
        Call /= L ;
        Put /= L ;
          
        //Averaging with blocking method
        A_Call += Call ;
        A2_Call += Call*Call ;
        sigma_Call = ( A2_Call / (double)(j+1.) ) - pow ( A_Call / (double)(j+1.), 2) ;
        if (j==0)
            sigma_Call = 0;
        else
            sigma_Call = std::sqrt (sigma_Call / (double)j ) ;
          
        A_Put += Put ;
        A2_Put += Put*Put ;
        sigma_Put = ( A2_Put / (double)(j+1.) ) - pow ( A_Put / (double)(j+1.), 2) ;
        if (j==0)
            sigma_Put = 0;
        else
            sigma_Put = std::sqrt (sigma_Put / (double)j ) ;
            
        //Writing on the file
        WriteDiscr << j << " " << A_Call/(j+1.) << " " << sigma_Call << " " << A_Put/(j+1.) << " " << sigma_Put << std::endl ;
            
    }
      
    WriteDiscr.close () ;
    
    
    gen.SaveSeed() ;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    return 0 ;
}
