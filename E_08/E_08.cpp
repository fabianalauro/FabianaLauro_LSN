//  LSN Exercise 5

//Optimisation of the ground state of a single quantum particle in a 1D space, through a Variational Monte Carlo technique.
//External potential: 𝑉(𝑥)=𝑥^4−5/2𝑥^2
//Approximate model of the wave function: Ψ𝜎,𝜇(𝑥) = 𝑒^(−(𝑥−𝜇)^2/2𝜎^2)+𝑒^(−(𝑥+𝜇)^2/2𝜎^2)

#include <iostream> //cout
#include <fstream> //ofstream
#include <cmath> //tanh

#include "Probability.hpp"
#include "Metropolis.hpp"
#include "random.h"
#include "BlockingAve.hpp"


double ene (double x, double mu, double sigma);

int main ()
{
    //Simulation parameters
    double mu = 0.836 ;
    double sigma = 0.624 ;
    const int wd = 20 ;
    double steplength = 1.2 ;
    
    //Data blocking parameters
    unsigned int N = 1E2 ;
    unsigned int M = 1E5;
    
    //Declaring my classes
    Random gen("Primes", "seed.in") ;
    Random * p_gen = &gen ;
       
    Variational Prob (sigma, mu) ;
    Basic_PDF * p_Prob = &Prob ;
       
    Metropolis1D METR ( 10., p_gen, p_Prob ) ;
    
    Blocking AVE (M, N) ;
    
    //Printing description
    std::cout << "Optimisation of the ground state of a single quantum particle in a 1D space" << std::endl;
    std::cout << "Using a Variational Monte Carlo technique" << std::endl << std::endl;
    
    std::cout << "External potential: V(x) = x^4 − 5/2 x^2" << std::endl;
    std::cout << "Approximate model of the wave function: Ψ𝜎,𝜇(x) = e^[−(𝑥−𝜇)^2/(2𝜎^2)]+e^[−(𝑥+𝜇)^2/(2𝜎^2)]" << std::endl << std::endl;
    
    std::cout << "ℏ = 1, m = 1" << std::endl << std::endl ;
    
    //Equilibration
    std::ofstream Out ("Equilibration.dat");
    std::cout << "EQUILIBRATION" << std::endl ;
    double alpha = 0.0 ;
    
    for (int i = 0; i < 100; i++)
    {
        METR.MRT2_Unif(steplength) ;
        Out << i+1 << std::setw(wd) << METR.x() << std::endl ;
        alpha += METR.alpha() ;
        
    }
    alpha /= 150. ;
    std::cout << "<alpha>=" << alpha << std::endl << std::endl ;
    
    Out.close() ;
    
    //Rough research
    std::cout << "Rough research" << std::endl;
    Out.open ("Lattice_rough.dat") ;
    steplength = 2.5 ;
    double best[3] = { AVE.ave(), sigma, mu };
    int sign = 1 ;
    mu = 0.05 ;
    sigma = 0.05 ;
    double step = 0.05 ;
    double sigma_max = 2. ;
    double mu_min = 0.05 ;
    double mu_max = 3. ;
    int counter = 0 ;

    while (sigma <= sigma_max)
    {
        while (mu >= mu_min && mu <= mu_max)
        {
            AVE.Reset() ;
            Prob.set_parameters (sigma, mu);
            for (int i = 0; i < M; i++)
            {
                METR.MRT2_Unif (steplength) ;
                double x = METR.x() ;
                alpha += METR.alpha() ;
                counter ++ ;
                double E = ene (x, mu, sigma);
                AVE.Add (E) ;
            }
            if (AVE.ave() < best[0])
            {
                best [0] = AVE.ave() ;
                best [1] = sigma;
                best [2] = mu ;
            }
            Out << sigma << std::setw(wd) << mu << std::setw(wd) << AVE.ave() << std::setw(wd) << AVE.err() << std::endl ;
            mu += sign*step ;
        }
        sign *= -1 ;
        mu += sign*step ;
        sigma += step ;
    }
    
    alpha /= counter ;
    std::cout << "<alpha>: " << alpha << std::endl ;
    std::cout << "BEST VALUES:" << std::endl;
    std::cout << "E=" << best[0] << " sigma=" << best[1] << " mu=" << best[2] << std::endl << std::endl ;
    Out.close() ;

   //Fine research
    std::cout << "Fine research" << std::endl;
    Out.open ("Lattice_fine.dat") ;
    steplength = 2.7 ;
    sign = 1 ;
    mu = best[2] - 0.05 ;
    sigma = best[1] - 0.05 ;
    sigma_max = best[1] + 0.05 ;
    mu_min = best[2] - 0.05 ;
    mu_max = best[2] + 0.05 ;
    step = 0.002 ;
    alpha = 0. ;
    counter = 0;
    
    std::ofstream Pos ("Positions.dat");
    while (sigma <= sigma_max)
    {
        while (mu >= mu_min && mu <= mu_max)
        {
            AVE.Reset() ;
            Prob.set_parameters (sigma, mu);
            for (int i = 0; i < M; i++)
            {
                METR.MRT2_Unif (steplength) ;
                double x = METR.x() ;
                if (sigma < 0.610001 && sigma > 0.60999 && mu < 0.82401 && mu > 0.82399)
                {
                    Pos << x << std::endl ;
                }
                alpha += METR.alpha() ;
                counter ++ ;
                double E = ene (x, mu, sigma);
                AVE.Add (E) ;
            }
            if (AVE.ave() < best[0])
            {
                best [0] = AVE.ave() ;
                best [1] = sigma;
                best [2] = mu ;
            }
            Out << sigma << std::setw(wd) << mu << std::setw(wd) << AVE.ave() << std::setw(wd) << AVE.err() << std::endl ;
            mu += sign*step ;
        }
        sign *= -1 ;
        mu += sign*step ;
        sigma += step ;
    }
    alpha /= counter ;
    std::cout << "<alpha>: " << alpha << std::endl ;
    std::cout << "BEST VALUES:" << std::endl;
    std::cout << "E=" << best[0] << " sigma=" << best[1] << " mu=" << best[2] << std::endl << std::endl ;
    
    Out.close() ;
    Pos.close() ;
    
    //<H>
    AVE.Reset() ;
    std::cout << "Calculation of <H> for the ground state" << std::endl ;
    Out.open("Average.dat") ;
    mu = best[2] ;
    sigma = best[1];
    alpha = 0.0 ;
    steplength = 2.8 ;
    for (int i = 0; i < M; i++)
    {
        METR.MRT2_Unif(steplength) ;
        double x = METR.x() ;
        alpha += METR.alpha() ;
        double E = ene (x, mu, sigma); 
        AVE.AddProg (E, Out) ;
    }
    alpha /= M ;
    std::cout << "<alpha>: " << alpha << std::endl << std::endl ;
       
    Out.close();
  
    return 0;
}




double ene (double x, double mu, double sigma)
   {
       return ( ( -1. / (2. * pow (sigma, 4)) ) * ( x*x + mu*mu - sigma*sigma - 2. * mu * x * tanh(mu*x/(sigma*sigma)) ) ) + pow(x, 4) - 5./2. * pow(x, 2);
   }
