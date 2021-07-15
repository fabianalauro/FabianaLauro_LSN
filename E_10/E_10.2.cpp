#include <iostream>
#include <vector>
#include <cmath>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include "random.h"
#include "Objects.hpp"
#include "mpi.h"


int main (int argc, char ** argv)
{
    //Parameters
    int ncities = 32 ; //Total number of cities in one trip
    int root = 0 ; //Rank of the root core
    int ntrips = 100 ; //Total number of trips in the population
    int ngens = 2000 ; //Number of generations
    int nmigr = 200; //Migrating between continents every nmigr generations

    //Parallelising
    int size, rank;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == root)
    {
        std::cout << "Program to solve the Travelling Salesman problem with a genetic algorithm, using parallel computing" << std::endl;
        std::cout << "Working with " << ncities << " cities " ;
        std::cout << "randomly placed inside a square" << std::endl;
        std::cout << "Population composed by " << ntrips << " trips" << std::endl;
        std::cout << "Total number of generations: " << ngens << std::endl << std::endl ;

        std::cout << "Using SPMD model " << std::endl ;
        std::cout << "Number of cores: " << size << std::endl  ;
        std::cout << "Migraton between continents every " << nmigr << " generations " << std::endl << std::endl;
    }

    double p_c = 0.7; //Crossover probability
    double p_m = 0.07 ; //Mutation probability
    double p = 5 ; //Selection parameter

    //Objects
    Random RND ("Primes", "seed.in", rank);
    Random * p_RND = &RND ;

    Population POP (ntrips, p_RND, size) ;
    Trip AppoTrip (ncities, p_RND);    


    std::ofstream evo ;
    std::ofstream best ;
    std::ofstream ave ;
    
    evo.open ("P_Evolution" + std::to_string(rank) + ".out");
    best.open ("P_FirstPath" + std::to_string(rank) + ".out");
    ave.open ("P_Ave" + std::to_string(rank) + ".out");
    


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
     
    //Broadcasting initial positions
    double x_i [ncities];
    double y_i [ncities];
  
    if(rank == root)
    {
        //AppoTrip.Print();
        for(int i=0; i<ncities; i++)
        {
            x_i[i] =  AppoTrip.cities()[i].x() ;
            y_i[i] =  AppoTrip.cities()[i].y() ;

        }
    }

    MPI_Bcast(x_i, ncities, MPI_REAL8, 0, MPI_COMM_WORLD);
    MPI_Bcast(y_i, ncities, MPI_REAL8, 0, MPI_COMM_WORLD);

    if(rank != root)
    {
        for(int i=0; i<ncities; i++)
        {
            AppoTrip.SetPositions (x_i, y_i) ;
        }
    }


    AppoTrip.measure() ;
    POP.AddTrip(AppoTrip, 0);
    POP.SetBest(AppoTrip);

    if(rank==0) POP.PrintSol (best) ;
    if(rank==1) POP.PrintSol (best) ;
    if(rank==2) POP.PrintSol (best) ;
    if(rank==3) POP.PrintSol (best) ;
    
    MPI_Barrier(MPI_COMM_WORLD);

    best.close();
    best.open ("P_BestPath" + std::to_string(rank) + ".out");

    
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
            std::cerr << rank << std::endl;   
            AppoTrip.Print();
            break;
        }
    }


    int migr_order [size];
    for (int i=0; i<size; i++)
        migr_order[i] = POP.migr()[i] ;
    

    //Evolution with migration
     for (int i=0; i<ngens; i++)
    {
        POP.Sort (true) ;
        POP.Selection(p);
        POP.GiveBirth(p_c);
        POP.Mutation(p_m);
        POP.NewGen();
        POP.Sort(true);
        //POP.PrintBest();
        POP.PrintBest(evo) ;
        POP.PrintAve (ave) ;
        if ((i+1)%nmigr == 0)
        {
            if (rank == root)
                std::cout << "MIGRATION #" << (i+1)/nmigr << std::endl ;
            if (rank == root)
            {
                for (int i=0;i<10;i++)
                    POP.SwapMigr();
                for (int i=0; i<size; i++)
                    migr_order [i] = POP.migr()[i] ;
            }
            MPI_Bcast(migr_order, 4 ,MPI_INTEGER, 0 , MPI_COMM_WORLD);
        
            POP.SetMigr (migr_order) ;
            

            MPI_Barrier(MPI_COMM_WORLD);

            POP.Migration(rank) ;

            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == root)
                std::cout <<"Done!" << std::endl << std::endl ;


        }   
    }





    POP.PrintSol(best);



    MPI_Finalize();

    evo.close();
    best.close();
    ave.close() ;

    p_RND->SaveSeed() ; 
    return 0;
}