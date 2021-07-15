

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
    if (argc != 4)
    {
        std::cerr << "Usage: <" << argv[0] << "> <circ> <p_c> <p_m> " << std::endl;
        return 1;
    }

    std::cout << "Program to solve the Travelling Salesman problem using a genetic algorithm" << std::endl;
    //Parameters
    bool circ = atoi(argv[1]); //if 1 works with 32 points on a circumference, if 0 works with 32 points inside a square
    int ncities = 32 ; //Total number of cities in one trip
    int ntrips = 100 ; //Total number of trips in the population
    int ngens ; //Number of generations
    if (circ == true)
        ngens = 500 ;
    else
        ngens = 5000;
    std::cout << "Working with " << ncities << " cities " ;
    if (circ == true)
        std::cout << "randomly placed on a circumference" << std::endl;
    else
        std::cout << "randomly placed inside a square" << std::endl;
    std::cout << "Population composed by " << ntrips << " trips" << std::endl;
    std::cout << "Total number of generations: " << ngens << std::endl << std::endl ;

    double p_c = atof(argv[2]); //Crossover probability
    double p_m = atof(argv[3]) ; //Mutation probability
    double p = 5 ; //Selection parameter

    //Objects
    Random RND ("Primes", "seed.in");
    Random * p_RND = &RND ;
    
    Population POP (ntrips, p_RND) ;
    Trip AppoTrip (ncities, p_RND);    

    std::ofstream evo ;
    std::ofstream best ;
    std::ofstream ave ;
    if (circ == true)
    {
        evo.open ("Circ_Evolution_" + std::to_string(p_c) + "_" + std::to_string(p_m) + ".out");
        best.open ("Circ_FirstPath_" + std::to_string(p_c) + "_" + std::to_string(p_m) + ".out");
        ave.open ("Circ_Ave_" + std::to_string(p_c) + "_" + std::to_string(p_m) + ".out");
    }
    else 
    {
        evo.open ("Square_Evolution_" + std::to_string(p_c) + "_" + std::to_string(p_m) + ".out");
        best.open ("Square_FirstPath_" + std::to_string(p_c) + "_" + std::to_string(p_m) + ".out");
        ave.open ("Square_Ave_" + std::to_string(p_c) + "_" + std::to_string(p_m) + ".out");
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
    POP.AddTrip(AppoTrip, 0);

    //Printing first configuration
    POP.PrintSol (best) ;
    best.close();

    if (circ == true)
        best.open ("Circ_BestPath_" + std::to_string(p_c) + "_" + std::to_string(p_m) + ".out");
    else
        best.open ("Square_BestPath_" + std::to_string(p_c) + "_" + std::to_string(p_m) + ".out");

    //Filling the population
    for (int i=1; i<ntrips; i++)
    {
        for (int j=0; j<50; j++)
           AppoTrip.Swap() ;
        if (AppoTrip.CheckFirst()==true)  
        {   
            AppoTrip.measure() ;
            POP.AddTrip(AppoTrip, i);
        }
        else
        {
            std::cerr << "ERROR: you have moved the first city. >:(" << std::endl ;    
            break;
        }
    }

    POP.SetBest(AppoTrip);
    //Evolving my population
    for (int i=0; i<ngens; i++)
    {
        POP.Sort (true) ;
        POP.Selection(p);
        POP.GiveBirth(p_c);
        POP.Mutation(p_m);
        POP.NewGen();
        POP.Sort(true);
        POP.PrintBest();
        POP.PrintBest(evo) ;
        POP.PrintAve (ave) ;
    }
    POP.PrintSol(best);
    //POP.PrintLength (true) ;
    //POP.Print(5);

    evo.close();
    best.close();
    ave.close() ;

    p_RND->SaveSeed() ; 
    return 0;
}

