//LSN exercise 02.2
//Simulating 3D random Walks on a cubic lattice and in the continuum

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"
#include "Position.hpp"

int main ()
{
    
    const unsigned int M = 1E4 ; //Number of throws
    const unsigned int N = 1E2 ; //Number of blocks
    unsigned int L = M / N ; //Number of throws in each block
    const unsigned int nstep = 100 ;  //Total number of the steps
    const double a = 1. ;   //Size of the step
    
    std::vector<double> Dr2 (nstep, 0.) ; //(DISCRETE) Sum of the r^2 inside a single block, for every value of the step
    std::vector<double> DA (nstep, 0.) ;  //(DISCRETE) Sum of the r^2 averaged on every block, for every value of the step
    std::vector<double> DA2 (nstep, 0.) ;  //(DISCRETE) Sum of the square of (r^2), with r^2 averaged on every block, for every value of the step
    std::vector<double> Dsigma  (nstep, 0.) ; //(DISCRETE) Standard deviation
    
    std::vector<double> Cr2 (nstep, 0.) ; //(CONTINUUM) Sum of the r^2 inside a single block, for every value of the step
       std::vector<double> CA (nstep, 0.) ;  //(CONTINUUM) Sum of the r^2 averaged on every block, for every value of the step
       std::vector<double> CA2 (nstep, 0.) ;  //(CONTINUUM) Sum of the square of (r^2) with r^2 averaged on every block, for every value of the step
       std::vector<double> Csigma  (nstep, 0.) ; //(CONTINUUM) Standard deviation
       
    //Constructing my Position classes, one for DISCRETE case and one for CONTINUUM case
    Position Dpos ;
    Position Cpos ;
    
    Random gen ("Primes", "seed.in") ;   //Random numbers generator
    
    //(DISCRETE) Variables to choose the directon
    double axis = 0. ;
    double direction = 0. ;
    
    //(CONTINUUM) Variables to choose the directon
    double teta = 0. ;
    double phi = 0. ;
    
    for (int j = 0; j < N; j ++)
    {
        for ( int k = 0; k < L; k ++ )
        {
            for ( int i = 0; i < nstep; i ++ )
            {
                //(DISCRETE) Making a step
                axis = gen.Rannyu (0., 3.) ;
                direction = gen.Rannyu () ;
        
                if (axis < 1.)
                {
                    if (direction < 0.5)
                        Dpos.StepCart (-a, 0., 0.) ;
                    else
                        Dpos.StepCart (a, 0., 0.) ;
                }
                else if (axis >= 2.)
                {
                    if (direction < 0.5)
                        Dpos.StepCart (0., 0., -a) ;
                    else
                        Dpos.StepCart (0., 0., a) ;
                }
                else
                {
                    if (direction < 0.5)
                        Dpos.StepCart (0., -a, 0.) ;
                    else
                        Dpos.StepCart (0., a, 0.) ;
                }
        
                Dr2[i] += Dpos.r2()  ;
                
                //(CONTINUUM) Making a step
                teta = acos (1 - 2 * gen.Rannyu ()) ;
                phi = 2 * M_PI * gen.Rannyu () ;
        
                Cpos.StepSph(phi, teta, a);
                
                Cr2[i] += Cpos.r2()  ;
                
            }
            //(DISCRETE) Re-initialising position
            Dpos.Set_x (0.) ;
            Dpos.Set_y (0.) ;
            Dpos.Set_z (0.) ;
            
            //(CONTINUUM) Re-initialising position
            Cpos.Set_x (0.) ;
            Cpos.Set_y (0.) ;
            Cpos.Set_z (0.) ;
        
        }
        //(DISCRETE) Averaging on blocks
        for ( int i = 0; i < nstep; i ++ )
        {
            DA[i] += Dr2 [i] / (double)L ;
            DA2[i] += pow ( Dr2 [i] / (double)L, 2 ) ;
            Dr2 [i] = 0. ;
        }
        
        //(CONTINUUM) Averaging on blocks
        for ( int i = 0; i < nstep; i ++ )
        {
            CA[i] += Cr2 [i] / (double)L ;
            CA2[i] += pow ( Cr2 [i] / (double)L, 2 ) ;
            Cr2 [i] = 0. ;
        }
        
    
    }
    //(DISCRETE) Sigma
    for ( int i = 0; i < nstep; i ++ )
    {
        Dsigma [i] = sqrt ( (DA2[i] / (double)N - ( pow (DA[i]/(double)N, 2) )) / (N-1.) ) ;
    }
    
    //(CONTINUUM) Sigma
    for ( int i = 0; i < nstep; i ++ )
    {
        Csigma [i] = sqrt ( (CA2[i] / (double)N - ( pow (CA[i]/(double)N, 2) )) / (N-1.) ) ;
    }
    
    //(DISCRETE) Writing the results on a file
    std::ofstream WriteDiscr ;
    WriteDiscr.open ( "RW_Discr.dat" ) ;
    
    if ( !WriteDiscr.is_open() )
    {
        std::cerr << "attention! Can't create output file for discrete random Walk" ;
        return 1 ;
    }
    
    for ( int i = 0; i < nstep; i ++ )
    {
        WriteDiscr << i << " " << sqrt (DA[i]/(double)N) << " " << (1. / ( 2. * sqrt (DA[i]) / N )) * Dsigma [i] << std::endl ;
    }
    
    WriteDiscr.close() ;
    
    //(CONTINUUM) Writing the results on a file
    std::ofstream WriteCont ;
      WriteCont.open ( "RW_Cont.dat" ) ;
      
      if ( !WriteCont.is_open() )
      {
          std::cerr << "attention! Can't create output file for continuous random Walk" ;
          return 1 ;
      }
      
      for ( int i = 0; i < nstep; i ++ )
      {
          WriteCont << i << " " << sqrt (CA[i]/(double)N) << " " << (1. / ( 2. * sqrt (CA[i]) / N )) * Csigma [i] << std::endl ;
      }
      
      WriteCont.close() ;
      
    
    return 0 ;
}
