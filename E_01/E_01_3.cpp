//LSN Exercise 02.1
//Simulating the Buffon's experimet


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"
#include "Fun_01.h"


//Make a picture of the estimation of  𝜋  and its uncertainty (Standard Deviation of the mean) with a large number of throws  𝑀  as a function of the number of blocks,  𝑁


int main ()
{
    
    //Distance between hoizontal lines, length of the needle
    const double d = 1. ;
    const double l = 0.8 ;
    
    //Number of throws, number of blocks, number of throws per block
    const int M = 1E6;
    const int N = 1E2 ;
    int L = M / N ;

    Random gen ("Primes", "seed.in") ;  //Random numbers generator
    
    double x_A = 0., x_B = 0., y_A = 0., y_B = 0., m = 0. ;  //Coordinates of the ends of the needle, slope of the needle
    int n_hit = 0, n_thr = 0 ;  //Counters for the throws and for the hits (intersection with a horiontal line is considered a hit)
    bool pos = true ;  //Variable to check on which side of A is the end B
    double pi = 0., pi2 = 0., sigma = 0. ; //Variables to estimate pi and its standard deviation
    
    //Creating the output file
    std::ofstream WritePi ;
    WritePi.open ("Pi.dat") ;
    if (!WritePi.is_open())
    {
        std::cerr << "Warning! Couldn't open Pi output file" ;
        return 1 ;
    }

    
    for (int j = 0; j < N; j++)
    {
        for (int i = 0; i < L; i++)
        {
            n_thr ++ ;
            
            //Generating position of end A of the needle in a 10x10 square
            x_A = gen.Rannyu (0., 10000.) ;
            y_A = gen.Rannyu (0., 10000.) ;
            
            //Generating direction of the needle in a circle of radius L, centered in A
            do
            {
                x_B = gen.Rannyu (x_A - l, x_A + l) ;
                y_B = gen.Rannyu (y_A - l, y_A + l) ;
            }
            while ( dist(x_A, y_A, x_B, y_B) > l ) ;
            
            //Calculating position of end B
            m = slope (x_A, y_A, x_B, y_B) ;
               
            if (y_B < y_A)
                pos = false ;
               
            y_B = y_cf (m, l, pos) + y_A ;
            
            //Checking if the needle intersects a horiziontal line
            if ( std::trunc(y_A) != std::trunc (y_B) )
                n_hit ++ ;
               
            }
        
        //Evaluating pi and writing the results on the file
        pi += (2 * l * n_thr) / (n_hit * d) ;
        pi2 += pow ((2 * l * n_thr) / (n_hit * d), 2) ;
        
        sigma = ( pi2 / (double)(j+1) ) - ( pi / (double)(j+1)*( pi / (double)(j+1) ) ) ;
        if (j==0)
            {
                sigma = 0;
            }
            else
            {
                sigma = std::sqrt (sigma / (double)j ) ;
            }
        
        WritePi << j << " " << pi / (double)(j+1) << " " << sigma << std::endl ;
        
        n_thr = 0 ;
        n_hit = 0 ;
        
    }
    
    return 0 ;
}
