#include <iostream>
#include <vector>
#include <cmath>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "random.h"
#include "Objects.hpp"

int main (int argc, char ** argv)
{
   
    if (argc != 2)
    {
        std::cerr << "Usage: <" << argv[0] << "> <circ>" << std::endl;
        return 1;
    }

    std::cout << "Program to solve the Travelling Salesman problem using Simulated Annealing" << std::endl;
    
    //Parameters
    bool circ = atoi(argv[1]); //if 1 works with 32 points on a circumference, if 0 works with 32 points inside a square
    int ncities = 32 ; //Total number of cities in one trip
    double i_temp = 10 ; // initial temperature
    double f_temp = 0.01 ;
    double coeff = 0.99999; //Temperature coefficient


    std::cout << "Working with " << ncities << " cities " ;
    if (circ == true)
        std::cout << "randomly placed on a circumference" << std::endl;
    else
        std::cout << "randomly placed inside a square" << std::endl;
    std::cout << "Temperature goes from " << i_temp << " to " << f_temp << " multiplicating by " << coeff << std::endl;

    //Objects
    Random RND ("Primes", "seed.in");
    Random * p_RND = &RND ;
    
    Trip AppoTrip (ncities, p_RND);    

    std::ofstream evo ;
    std::ofstream evobest ;
    std::ofstream best ;
    if (circ == true)
    {
        evo.open ("Circ_Evolution.out", std::ios::app);
        evobest.open ("Circ_EvolutionBest.out", std::ios::app);
        best.open ("Circ_FirstPath.out");
    }
    else 
    {
        evo.open ("Square_Evolution.out", std::ios::app);
        evobest.open ("Square_EvolutionBest.out", std::ios::app);
        best.open ("Square_FirstPath.out");
    }

    if (circ == true)
    {
        //Creating ncities cities on a circumference with r=1
        for (int i=0; i<ncities; i++)
        {
            double theta = RND.Rannyu(0, 2*M_PI) ;
            if (i==0)
            {
                City C (i, true, cos(theta), sin(theta));
                AppoTrip.AddCity(C, i);
            }
            else
            {
                City C (i, false, cos(theta), sin(theta));
                AppoTrip.AddCity(C, i);
            } 
        }
    }
    else 
    {
        //Creating ncities cities inside a square with l=1
        for (int i=0; i<ncities; i++)
        {
            double x = RND.Rannyu(0, 1) ;
            double y = RND.Rannyu(0, 1) ;
            if (i==0)
            {
                City C (i, true, x, y);
                AppoTrip.AddCity(C, i);
            }
            else
            {
                City C (i, false, x, y);
                AppoTrip.AddCity(C, i);
            } 
        }
    }

    
    AppoTrip.measure() ;
    AppoTrip.Print (best) ;
    best.close();

    if (circ == true)
        best.open ("Circ_BestPath.out");
    else
        best.open ("Square_BestPath.out");


    Annealing SIM (ncities, AppoTrip, p_RND) ;
    double temp = i_temp ;
    double alpha; 
    //SIM.Best().Print () ;


    //for (int j = 0; j<200; j++)
    while (temp > f_temp)
    {   alpha = 0. ;
        SIM.SetBeta(1./temp) ;
        SIM.MRT2() ;
        std::cout << "Temperature=" << temp << "\t alpha=" << SIM.alpha() << "\t best length=" << SIM.Best().length() << std::endl;
        SIM.PrintEvo (evo) ; 
        SIM.PrintEvoBest (evobest);
        temp *= coeff ;
    }

    SIM.Print(best) ;

    evobest.close();
    evo.close () ;
    best.close();  

    p_RND->SaveSeed() ; 
    return 0;
}