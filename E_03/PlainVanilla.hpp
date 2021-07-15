//  Class to simulate a PlainVanilla (European) option in a financial market. Thanks to a random numbers generator, it can calculate the final asset price of the product using the Geometric Brownian Motion equation, and then calculate the profit for a Call or for a Put option.
//  It works on a statistic basis: make a lot of simulations to have good results!
//


#ifndef PlainVanilla_hpp
#define PlainVanilla_hpp

#include <cmath>
#include "random.h"

class PlainVanilla
{
private:
    Random * m_gen ; //Random numbers generator
    double m_t0, m_t1, m_K  ; //Option parameters: initial time, delivery time, strike price
    double m_S0, m_r, m_sigma ;  //Market parameters: initial asset price, risk-free interest rate, volatility
    double m_S1 ;  //Final asset price
    
public:
    
    PlainVanilla () ; //Constructor with no arguments: sets the parameters and m_S1 to 0 , the pointer to NULL
    PlainVanilla (double t0, double t1, double S0, double K, double r, double sigma, Random * gen) ; //Constructor that sets the parameters and the pointer to the given values
    
    void setRandom (Random * gen) ; //method to set m_gen
    void setParameters (double t0, double t1, double S0, double K, double r, double sigma) ;  //Method to set the parameters
    void setS1_dir () ;  //Method to calculate a new value of the final asset price in a direct way, and overwrite m_S1 with the new value
    void setS1_discr (unsigned int n_steps) ; //Method to calculate a new value of the final asset price through n_steps discrete steps, and overwrite m_S1 with the new value
    
    double S1 () ; //Returns the value of m_S1
    
    double Call_Profit () ;  //Calculates the profit for a Call option
    double Put_Profit () ;  //Calculates the profit for a Put option
    
    
} ;




#endif /* PlainVanilla_hpp */
