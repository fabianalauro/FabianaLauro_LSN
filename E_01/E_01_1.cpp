// LSN Exercise 01.1
// Testing the random numbers generator


#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"



int main ()
{
    
    const unsigned int M = 100000; //Number of throws
    const unsigned int N = 100;   //Number of blocks
    unsigned int L = (double) M / (double) N;   //Number of throws in each block
    
    Random gen ("Primes", "seed.in") ;  //random numbers generator
    

//1
//Making a picture of the estimation of  ⟨𝑟⟩  and its uncertainty with a large number of throws 𝑀, as a function of the number of blocks,  𝑁
    
    double m = 0 ;  //The random numbers generated
    double A = 0 ;  //Average in every block
    double A2 = 0 ;  //Square of the average in every block
    double sigma = 0 ;  //Standard Deviation
    
    //Creating output file
    std::ofstream WriteData;
    WriteData.open ("Data_01_1.dat") ;
    if (!WriteData)
    {
        std::cerr << "Warning! Couldn't open Average output file" ;
        return 1 ;
    }
    
    //Generating the random numbers, calculating the average with the blocks method and writing the results on the file
    for (int i=0; i < N; i++)
    {
        for (int j = 0; j < L; j++)
        {
            m += gen.Rannyu() ;
        }
        m /= L ;
        A += m ;
        A2 += m * m ;
        
        sigma = ( A2 / (double)(i+1) ) - ( A / (double)(i+1)*( A / (double)(i+1) ) ) ;
        
        if (i==0)
        {
            sigma = 0;
        }
        else
        {
            sigma = std::sqrt (sigma / i ) ;
        }
        
        WriteData << i << " " << (A / (double)(i+1)) << " " << sigma << std::endl ;
        m = 0 ;
    }
    
    WriteData.close() ;
    
    
//2
// Making a picture of the estimation of  𝜎^2  and its uncertainty with a large number of throws  𝑀, as a function of the number of blocks,  𝑁
    
    //Re-initialising
    A = 0 ;
    A2 = 0 ;
    sigma = 0 ;
    
    //Creating output file
    std::ofstream WriteSigma;
    WriteSigma.open ("Sigma_01_1.dat") ;
    if (!WriteSigma)
    {
        std::cerr << "Warning! Couldn't open Sigma output file" ;
        return 1 ;
    }
    
    //Generating the random numbers, calculating the standard deviation with the blocks method and writings the results on the file
    for (int i=0; i < N; i++)
       {
           for (int j = 0; j < L; j++)
           {
               m += std::pow (gen.Rannyu()-0.5, 2) ;
           }
           
           m /= L ;
           A += m ;
           A2 += m * m ;
           
           sigma = ( A2 / (double)(i+1) ) - ( A / (double)(i+1)*( A / (double)(i+1) ) ) ;
           
           if (i==0)
           {
               sigma = 0;
           }
           else
           {
               sigma = std::sqrt (sigma / i ) ;
           }
           
           WriteSigma << i << " " << A / (double)(i+1) << " " << sigma << std::endl ;
           m = 0 ;
       }
       
       WriteSigma.close() ;
    
//3
//Dividing  [0,1]  into  𝑀  identical sub-intervals and implementing the  𝜒2  test.
    
    
    const unsigned int W = 100;  //Number of sub-intervals (M has already been used)
    const unsigned int n = 10000 ; //Number of random numbers generated
    const unsigned int j = 100 ;  //Number of 𝜒2 tests I make
    double E = (double)n/(double)W ;
    
    std::vector<int> BINS (W, 0) ;  //Counter for every sub-interval
    double Chi2 = 0 ;  //Variable for Chi^2
    
    //Creating output file
    std::ofstream WriteChi ;
    WriteChi.open ("Chi_01_1.dat");
    if (!WriteChi)
    {
        std::cerr << "Warning! Couldn't open Chi2 output file" ;
        return 1 ;
    }
    
    //Generating the random numbers, calculating 𝜒2 and writing the results on the file
    for (int g = 0; g < j; g++)
    {
        for (int k = 0; k < n; k ++)
        {
            m = gen.Rannyu() * 100 ;
            m = trunc (m) ;
            BINS[m] ++ ;
        }
        for (int i = 0; i < W; i++)
        {
            Chi2 += (BINS[i] - E) * (BINS[i] - E) / E ;
        }
        WriteChi << g << " " << Chi2 << std::endl ;
        
        //Re-initialising
        Chi2 = 0 ;
        for (int i = 0; i < W; i++)
        {
            BINS[i] = 0 ;
        }
    }
    
    WriteChi.close() ;
    
    return 0 ;
}
