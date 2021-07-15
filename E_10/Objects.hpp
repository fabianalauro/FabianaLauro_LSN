//
//
//
//
//  Created by Fabiana on 06/07/21.
//

#ifndef Objects_hpp
#define Objects_hpp


#include <iostream>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "random.h"
#include "mpi.h"


//City
class City
{
    
private:
    
    int m_label ;
    bool m_first ;
    double m_x, m_y ;
    
protected:
    
public:
    
    City() ;
    City (int label, bool first, double x, double y) ;
    ~City() {};

    City& operator= (const City& C) ;
    
    //Return data members
    int label() const {return m_label;};
    bool first() const {return m_first;};
    double x() const {return m_x;};
    double y() const {return m_y;};

    //Set data members
    void SetX (double x);
    void SetY (double y);
    void SetL (int l) ;

    void PrintCity () ; //Prints on screen all the data members
    void PrintCity (std::ofstream &out) ; //Prints all the data members on a file

};


//Trip
class Trip
{
    
private:
    
    int m_ncities ; //Number of cities in this trip
    std::vector<City> m_cities ; //Vector containing the cities in the right order
    double m_length ; //length of the trip
    Random * m_rnd ; //Pointer to random number generator
    
protected:
    
public:
    
    Trip() ;
    Trip (int N, Random * gen ) ;
    //Trip (Trip& T) ;
    
    ~Trip() {};
    
    //Return data members
    int ncities () const {return m_ncities;};
    double length () const {return m_length;} ;
    Random * rnd() const {return m_rnd;};
    std::vector<City> cities() const {return m_cities;};
    
    //Set data members
    void SetCities (std::vector<City> V) ;
    void SetPositions (double * x, double * y) ;
    void SetPositions (int * l, double * x, double * y) ;
    
    Trip& operator= (const Trip& T) ;

    bool CheckFirst () ; //Checks if the first city in m_cities is the right one
    void AddCity (City C, int index) ; //Adds a city in m_cities, in the position indicated by the index given
    void measure ()  ; //Calculates the length of the trip, updates m_length
    int Pbc (int index, int size)  ; //m_cities[m_ncities] = m_cities[1] and so on (doesn't touch n_cities[0])
    void Print () ; //Calls PrintCity() for every city in m_cities
    void Print (std::ofstream &out) ; //Calls PrintCity() for every city in m_cities


    //Mutations
    void Swap () ; //Swaps two randomly chosen cities inside m_cities (doesn't touch the first)
    void SwapCont() ; //Swaps two contiguous cities inside m_cities (doesn't touch the first)
    void SwapM() ; //Swaps m contiguous cities with other m contiguous cities in m_cities (doesn't touch the first)
    void Shift() ; //Shifts forward m contiguous cities of n positions
    void Reverse() ; //Inverts of the order in which appear in the path m cities
    void Migration () ; //Performs migration between continents

};

//Annealing
class Annealing
{
    
    
private:
    Trip m_TheTrip ; //
    Trip m_NewTrip ; //
    Random * m_rnd ; //Pointer to random numbers generator
    double m_beta ; //beta=1/T
    int m_count ; //Progressive number of step
    int m_bestcount ; //Record of the step in which the best trip eas found
    double m_alpha ; //Records the acceptance rate
    Trip m_best ; //Records the best trip so far
    
protected:
        
public:
    
    Annealing () ;
    Annealing (int ncities, Trip T, Random * gen) ;
    ~Annealing () {};
    
    //Return data members 
    double alpha () {return m_alpha;};
    double beta () {return m_beta;};
    Trip Best() {return m_best;};

    //Utility
    void SetBeta (double beta) {m_beta = beta;};

    void Mutate () ;
    void MRT2 () ;

    void Print (bool The) ;

    void Print (std::ofstream &out) ;
    void PrintEvo (std::ofstream &out) ;
    void PrintEvoBest (std::ofstream &out) ;


};


//Population
class Population
{

private:
    int m_ntrips ; //Number of thrips in the population
    std::vector<Trip> m_trips ; //Vector containing the trips in the population
    std::vector<Trip> m_deliveryroom ; //Vector where the crossover happens
    Random * m_rnd ; //Pointer to random numbers generator
    int m_generation ; //Progressive number of generation
    int m_bestgen ; //Records the best generation so far
    Trip m_best ; //Records the best trip so far
    std::vector<int> m_migr  ; //Decides the migrations: cores at indexes 1,2 migrate to cores at indexes 3,4
    
protected:
        
public:
    
    Population () ;
    Population (int N, Random * gen) ;
    Population (int N, Random * gen, int ncores) ;
    
    //Return data members
    std::vector <Trip> trips() {return m_trips;};
    std::vector <int> migr () {return m_migr;};

    //Utility
    void SetBest (Trip T) {m_best=T;} ;
    void AddTrip (Trip T, int index) ; //Adds a trip to m_trips, in the index given
    void Swap (int index1, int index2, bool pop) ; //if pop=true swaps two trips in the vector m_trips, otherwise two trips in the vector m_deliveryroom
    void Sort (bool pop) ; //if pop=true sorts the trips in m_trips by their length, otherwise sorts the babies in m_deliveryroom by their length
    void SwapMigr () ; //Swaps two numbers in m_migr
    void SetMigr (int * v) ;

    //Evolution
    void Selection (double p); //Selects the parents of the new population
    void Crossover (Trip &P1, Trip &P2); //Implements half of the crossover (only for one child)
    void GiveBirth (double p) ; //Implements the crossover on all the elements in the deliveryroom
    void Mutation(double p) ; //Performs mutations on the delivery room, all with different probability
    void NewGen () ; //Creates a new generation keeping the best half of the old population and the best half of the new population
    void Migration (int rank); //Performs migration

    //Export Results
    void PrintLength (bool pop) ; //if pop=true prints the lengths of the elements in m_trips, otherwise prints the lengths of the elements in m_deliveryroom
    void Print () ;//Calls Print() for every trip in m_trips
    void Print (int n) ;//Calls Print() for the first n trips in m_trips
    void PrintBest(); //Prints m_generation, the best trip length and m_bestgen
    void PrintBest(std::ofstream &out); //Prints on file m_generation, the best trip length and m_bestgen
    void PrintAve (std::ofstream &out) ; //Prints on file m_generation and the average trip length of the best part of the population
    void PrintSol (std::ofstream &out); //Prints the solution (the cities in the right order) on a file
    void PrintSol (); //Prints on screen the solution (the cities in the right order) 

};


#endif /* Objects_hpp */

