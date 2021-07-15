//  LSN Exercise 5
//  
// Calculating the average radius of a Hydrogen atom, in its ground state and in one of the three  2𝑝  excited states


#include <iostream> //cin, cout, cerr: interact with standard input/output
#include <fstream>   //interact with files
#include <cmath>  //sqrt, pow
//#include <cstdlib>  //atoi

#include "random.h"
#include "Probability.hpp"
#include "Metropolis.hpp"
#include "BlockingAve.hpp"


int main (int argc, char ** argv)
{
    //(UNIFORM + GAUSS) Preparing and initialising Metropolis classes
    Random gen("Primes", "seed.in") ;
 
    Random * p_gen = &gen ;
    
    Hydr100 Wave1 ;
    Basic_PDF * p_Wave = &Wave1 ;
    
    Metropolis3D U_METR ( 20., 20., 20., p_gen, p_Wave ) ;
    Metropolis3D G_METR ( 20., 20., 20., p_gen, p_Wave ) ;
    

    //declaring constants for data blocking, creating Data Blocking class (UNIFORM)
    const unsigned int M = 1E6 ; //Number of throws
    const unsigned int N = 1E2 ;  //Number of blocks
    
    Blocking U_Ave (M, N) ;
    Blocking G_Ave (M, N) ;
    
    //Variables to check the rate of acceptance
    double Ualpha = 0. ;
    double Galpha = 0. ;
    
    //(UNIFORM + GAUSS) File streams
    std::ofstream WriteDataU ;
    std::ofstream WriteDataG ;
    std::ofstream WritePos ;
    //std::ofstream * p_WritePos = & WritePos ;
    
//1. psi_1,0,0
    
    //Equilibration: starting far from origin, getting closer (1000 steps)
    
    const unsigned int FFO_nsteps = 1000 ;
    
    WritePos.open ("FFO_100.dat") ;
    
    if ( !WritePos.is_open() )
    {
        std::cerr << "Attention! Can't create output file for 1,0,0 equilibration!" ;
        return 1 ;
    }
    
    for (int i = 0; i < FFO_nsteps; i++)
    {
        U_METR.MRT2_Unif (1.3) ;
        Ualpha += U_METR.alpha() ;
        
        G_METR.MRT2_Gauss (0.8) ;
        Galpha += G_METR.alpha() ;
        
        WritePos << U_METR.x() << " " << U_METR.y() << " " << U_METR.z() << " " << G_METR.x() << " " << G_METR.y() << " " << G_METR.z() << std::endl ;
        
    }
    
    Ualpha /= FFO_nsteps ;
    Galpha /= FFO_nsteps ;
    std::cout << "AVERAGE VALUES OF ALPHA DURING 1,0,0 EQUILIBRATION PATH:" << std::endl ;
    std::cout << Ualpha << " " << Galpha << std::endl << std::endl ;
    
    Ualpha = 0. ;
    Galpha = 0. ;
    
    WritePos.close() ;
    
    //(UNIFORM + GAUSS) Creating output files
    WriteDataU.open ( "r_100U.dat" ) ;
    WriteDataG.open ( "r_100G.dat" ) ;
    WritePos.open ( "pos_100.dat" ) ;
       
    if ( !WriteDataU.is_open() )
    {
        std::cerr << "Attention! Can't create output file for UNIFORM 1,0,0 results!" ;
        return 1 ;
    }
    
    if ( !WriteDataG.is_open() )
    {
        std::cerr << "Attention! Can't create output file for GAUSS 1,0,0 results!" ;
        return 1 ;
    }
    
    if ( !WritePos.is_open() )
    {
        std::cerr << "Attention! Can't create output file for 1,0,0 positions!" ;
        return 1 ;
    }
       
    //(UNIFORM + GAUSS) Calculating <r>, writing positions every 100 steps
    
    for ( int j = 0; j < M; j ++ )
    {
        U_METR.MRT2_Unif (1.2) ;
        U_Ave.Add ( U_METR.r(), WriteDataU ) ;
        Ualpha += U_METR.alpha() ;
        
        G_METR.MRT2_Gauss (0.8) ;
        G_Ave.Add ( G_METR.r(), WriteDataG ) ;
        Galpha += G_METR.alpha() ;
        
        if (j % 100 == 0)
            WritePos << U_METR.x() << " " << U_METR.y() << " " << U_METR.z() << " " << G_METR.x() << " " << G_METR.y() << " " << G_METR.z() << std::endl ;
        
    }
  
    WriteDataU.close() ;
    WriteDataG.close() ;
    WritePos.close() ;
    
    Ualpha /= M ;
    Galpha /= M ;
    std::cout << "AVERAGE VALUES OF ALPHA DURING 1,0,0 SIMULATION PATH:" << std::endl ;
    std::cout << Ualpha << " " << Galpha << std::endl << std::endl  ;

/*------------------------------------------------------------------------------------*/
    
//2. psi_2,1,0
    
    //Re-initialising
    
    Hydr210 Wave2 ;
    p_Wave = &Wave2 ;
    
    U_METR.Set_coord (20., 20., 20.) ;
    U_METR.Set_PDF (p_Wave) ;
    
    G_METR.Set_coord (20., 20., 20.) ;
    G_METR.Set_PDF (p_Wave) ;
    
    U_Ave.Reset() ;
    G_Ave.Reset() ;
    
    std::cout << "--------------------------------------------------------------------------" << std::endl << std::endl ;

    //Equilibration: starting far from origin, getting closer (1000 steps)
    
    WritePos.open ("FFO_210.dat") ;
    
    if ( !WritePos.is_open() )
    {
        std::cerr << "Attention! Can't create output file for 2,1,0 equilibration path!" ;
        return 1 ;
    }
    
    for (int i = 0; i < FFO_nsteps; i++)
    {
        U_METR.MRT2_Unif (2.9) ;
        Ualpha += U_METR.alpha() ;
        
        G_METR.MRT2_Gauss (1.9) ;
        Galpha += G_METR.alpha() ;
        
        WritePos << U_METR.x() << " " << U_METR.y() << " " << U_METR.z() << " " << G_METR.x() << " " << G_METR.y() << " " << G_METR.z() << std::endl ;
        
    }
    
    Ualpha /= FFO_nsteps ;
    Galpha /= FFO_nsteps ;
    std::cout << "AVERAGE VALUES OF ALPHA DURING THE 2,1,0 EQUILIBRATION PATH:" << std::endl ;
    std::cout << Ualpha << " " << Galpha << std::endl << std::endl ;
    
    Ualpha = 0. ;
    Galpha = 0. ;
    
    WritePos.close() ;
    
    
    //(UNIFORM + GAUSS) Creating output files
    WriteDataU.open ( "r_210U.dat" ) ;
    WriteDataG.open ( "r_210G.dat" ) ;
    WritePos.open ( "pos_210.dat" ) ;
       
    if ( !WriteDataU.is_open() )
    {
        std::cerr << "Attention! Can't create output file for UNIFORM 2,1,0 wave function results!" ;
        return 1 ;
    }
    
    if ( !WriteDataU.is_open() )
    {
        std::cerr << "Attention! Can't create output file for GAUSS 2,1,0 wave function results!" ;
        return 1 ;
    }
       
    if ( !WritePos.is_open() )
    {
        std::cerr << "Attention! Can't create output file for 2,1,0 wave function positions!" ;
        return 1 ;
    }
    
    for ( int j = 0; j < M; j ++ )
    {
        U_METR.MRT2_Unif (2.8) ;
        U_Ave.Add ( U_METR.r(), WriteDataU ) ;
        Ualpha += U_METR.alpha() ;
        
        G_METR.MRT2_Gauss (1.9) ;
        G_Ave.Add ( G_METR.r(), WriteDataG ) ;
        Galpha += G_METR.alpha() ;
        
        if (j % 100 == 0)
            WritePos << U_METR.x() << " " << U_METR.y() << " " << U_METR.z() << " " << G_METR.x() << " " << G_METR.y() << " " << G_METR.z() << std::endl ;
        
    }
         
    WriteDataU.close() ;
    WriteDataG.close() ;
    WritePos.close() ;
    
    Ualpha /= M ;
    Galpha /= M ;
    std::cout << "AVERAGE VALUES OF ALPHA DURING THE 2,1,0 SIMULATION PATH:" << std::endl ;
    std::cout << Ualpha << " " << Galpha << std::endl << std::endl  ;
    

    
    return 0 ;
}



