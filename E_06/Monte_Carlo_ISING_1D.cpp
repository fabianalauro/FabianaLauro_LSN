/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Simulation of a 1D Ising model, using Metropolis and Gibbs algorithms

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(int argc, char ** argv)
{
    if (argc != 2 )
    {
        cerr << "Usage: <" << argv[0] << "> <num> (temp = 0.5 + (0.15*num))" << endl ;
        exit (-1) ;
    }
    num = atoi(argv[1]) ;
    Input(); //Initialisation
    
    for (int i = 0; i < 200; i++) //Equilibration
    {
        Move (metro) ;
        Measure () ;
        if (temp==0.5 || temp==2. || temp==1.25)
            PrintEqui () ;
    }
    //Saving config at the end of the equilibration (only at some temperatures)
    if (temp==0.5 || temp==2. || temp==1.25)
    {
        ofstream conf ;
        if (metro==1)
            conf.open ("Metro/config.f_T" + to_string(temp) + "_M" + to_string(h)) ;
        else
            conf.open ("Gibbs/config.f_T" + to_string(temp) + "_M" + to_string(h)) ;
        for (int i=0; i<nspin; ++i)
        {
            conf << s[i] << endl;
        }
        conf.close();
    }
    
    for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
    {
        Reset(iblk);   //Reset block averages
        for(int istep=1; istep <= nstep; ++istep)
        {
            Move(metro);
            Measure();
            Accumulate(); //Update block averages
        }
        Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
    ifstream ReadInput;

    cout << "Classic 1D Ising model             " << endl;
    cout << "Monte Carlo simulation             " << endl << endl;
    cout << "Nearest neighbour interaction      " << endl << endl;
    cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
    cout << "The program uses k_B=1 and mu_B=1 units " << endl;

    //Read seed for random numbers
    int p1, p2;
    ifstream Primes("Primes");
    Primes >> p1 >> p2 ;
    Primes.close();

    ifstream input("seed.out");
    input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    rnd.SetRandom(seed,p1,p2);
    input.close();
  
    //Read input informations
    ReadInput.open("input.dat");

    //ReadInput >> temp;
    temp = 0.5 + (0.15*num) ;
    beta = 1.0/temp;
    cout << "Temperature = " << temp << endl;

    ReadInput >> nspin;
    cout << "Number of spins = " << nspin << endl;

    ReadInput >> J;
    cout << "Exchange interaction = " << J << endl;

    ReadInput >> h;
    cout << "External field = " << h << endl << endl;
    
    ReadInput >> metro; // if = 1 Metropolis else Gibbs
    
    ReadInput >> restart;

    ReadInput >> nblk;

    ReadInput >> nstep; //Number of steps inside one block

    if(metro == 1) cout << "The program performs Metropolis moves" << endl;
    else cout << "The program performs Gibbs moves" << endl;
    
    cout << "Calculating average values with the blocking method" << endl ;
    cout << "Number of blocks = " << nblk << endl;
    cout << "Number of steps in one block = " << nstep << endl << endl;
    ReadInput.close();


    //Prepare arrays for measurements
    iu = 0; //Energy
    ic = 1; //Heat capacity
    im = 2; //Magnetisation
    ix = 3; //Magnetic susceptibility
    
    n_props = 4; //Number of observables

    //initial configuration (set randomly)
    if (restart == 0)
    {
        cout << "Setting the initial configuration randomly (T = infty)" << endl << endl ;
        for (int i=0; i<nspin; ++i)
        {
            if(rnd.Rannyu() >= 0.5)
                s[i] = 1;
            else
                s[i] = -1;
        }
    }
    else
    {
        cout << "Importing previous configuration from file config.final" << endl << endl ;
        ifstream conf;
        conf.open ("config.final") ;
        if (conf.fail())
        {
            cerr << "Problems opening the file config.out!" << endl ;
            exit (1) ;
        }
        for (int i=0; i<nspin; ++i)
        {
            conf >> s[i] ;
        }
        conf.close() ;
    }
    
    //Saving initial configuration (only at some temperatures)
    if (temp==0.5 || temp==2. || temp==1.25)
    {
        ofstream conf ;
        if (metro==1)
            conf.open ("Metro/config.i_T" + to_string(temp) + "_M" + to_string(h)) ;
        else
            conf.open ("Gibbs/config.i_T" + to_string(temp) + "_M" + to_string(h)) ;
        for (int i=0; i<nspin; ++i)
        {
            conf << s[i] << endl;
        }
        conf.close();
        
    }
  
    //Evaluate energy etc. of the initial configuration
    Measure();

    //Print initial values for the potential energy and virial
    cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro)
{
    int o;

    for (int i=0; i<nspin; ++i)
    {
        //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
        o = (int)(rnd.Rannyu()*nspin);

        if (metro == 1) //Metropolis
        {
            double alpha = fmin ( 1, Probability (beta, Boltzmann (s[o], o) - Boltzmann (-s[o], o)) ) ;
            double r = rnd.Rannyu() ;
            attempted ++ ;
            if (r < alpha || r == alpha)
            {
                s[o] *= -1 ;
                accepted ++ ;
            }
        }
        else //Gibbs sampling
        {
            double prob = 1. /( 1. + Probability ( beta, ( Boltzmann(1, o) - Boltzmann(-1, o) ) ) );
            double r = rnd.Rannyu() ;
            if (r < prob || r == prob)
                s[o] = 1 ;
            else
                s[o] = -1 ;
        }
    }
}

double Boltzmann(int sm, int ip)
{
    double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
    return ene;
}

double Probability (double beta, double E)
{
    return exp (beta * E) ;
}


void Measure()
{
    //int bin;
    double u = 0.0, m = 0.0;

//cycle over spins
    for (int i=0; i<nspin; ++i)
    {
        u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);  // sommo i valori di H per poi calcolare <H>
        m += s[i] ; //Sommo i valori di s
    }
    walker[iu] = u ; //Sto sommando gli H
    walker[ic] = u*u ; //Sto sommando gli H^2
    walker[im] = m ; //Sto sommando gli s
    walker [ix] = m * m ; //Sto sommando gli s^2
}

void PrintEqui (void)
{
    ofstream Equi ;
    if (metro==1)
        Equi.open ("Metro/equi.ene_T" + to_string(temp) + "_M" + to_string(h), ios::app)  ;
    else
        Equi.open ("Gibbs/equi.ene_T" + to_string(temp) + "_M" + to_string(h), ios::app)  ;
    Equi << walker[iu]/(double)nspin << endl ;
    Equi.close();
}


void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) //Update block averages
{

    for(int i=0; i<n_props; ++i)
    {
       blk_av[i] = blk_av[i] + walker[i];
    }
    blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;
   const int wd = 20;
    
    cout << "Block number " << iblk << endl;
    if (metro == 1)
        cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    //Internal energy: averaging H/N
    if (metro==1)
        Ene.open("Metro/output.ene_M" + to_string(h),ios::app);
    else
        Ene.open("Gibbs/output.ene_M" + to_string(h),ios::app);
    stima_u = blk_av[iu] / blk_norm; //<H> inside the block
    glob_av[iu]  += stima_u; //A for <H>
    glob_av2[iu] += stima_u * stima_u; //A2 for <H>
    err_u = Error (glob_av[iu], glob_av2[iu], iblk) ;  //Sigma for <H>
    
    Ene << setw(wd) << iblk <<  setw(wd)
        << stima_u/(double)nspin << setw(wd) << glob_av[iu]/(double)iblk/(double)nspin << setw(wd)
        << err_u/(double)nspin << endl;
    Ene.close();
    
    //Heat capacity: averaging H^2
    if (metro==1)
           Heat.open("Metro/output.heat_M" + to_string(h),ios::app);
       else
           Heat.open("Gibbs/output.heat_M" + to_string(h),ios::app);    stima_c = blk_av[ic] / blk_norm ; // <H^2> inside the block
    glob_av[ic] += stima_c ; //A for <H^2>
    glob_av2[ic] += stima_c * stima_c ; //A2 for <H^2>
    err_c = Error (glob_av[ic], glob_av2[ic], iblk) ; //Sigma for <H^2>
    
    Heat << setw(wd) << iblk << setw(wd)
         << beta*beta * (stima_c - pow (stima_u, 2)) << setw(wd)
         << beta*beta * (glob_av[ic]/(double)iblk - pow (glob_av[iu]/(double)iblk, 2))/(double)nspin << setw(wd)
         << beta * beta * sqrt ( pow ( err_c, 2 ) + pow ( 2 * fabs(glob_av[iu])/(double)iblk * err_u, 2 ) )/(double)nspin << endl ;
    Heat.close() ;
    
    //Magnetisation: averaging s
    if (metro==1)
           Mag.open("Metro/output.mag_M" + to_string(h),ios::app);
       else
           Mag.open("Gibbs/output.mag_M" + to_string(h),ios::app);    stima_m = blk_av[im] / blk_norm ; //<s> inside the block
    glob_av[im] += stima_m ; //A for <s>
    glob_av2[im] += stima_m * stima_m ; //A2 for <s>
    err_m = Error (glob_av[im], glob_av2[im], iblk) ; //Sigma for <s>
    
    Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd)
        << glob_av[im]/(double)iblk/(double)nspin << setw(wd)
        << err_m/(double)nspin << endl;
    Mag.close();
    
    //Magnetic susceptibility: averaging s^2
    if (metro==1)
           Chi.open("Metro/output.chi_M" + to_string(h),ios::app);
       else
           Chi.open("Gibbs/output.chi_M" + to_string(h),ios::app);    stima_x = blk_av[ix] / blk_norm ; //<s^2> inside the block
    glob_av[ix] += stima_x ; //A for <s^2>
    glob_av2[ix] += stima_x*stima_x ; //A2 for <s^2>
    err_x = Error (glob_av[ix], glob_av2[ix], iblk) ; //Sigma for <s^2>
    
    if (h==0)
        Chi << setw(wd) << iblk << setw(wd) << beta * (stima_x) << setw(wd)
            << beta * glob_av[ix]/(double)iblk/(double)nspin  << setw(wd)
            << beta * err_x /(double)nspin << endl ;
    else
        Chi << setw(wd) << iblk << setw(wd) << beta * (stima_x) << setw(wd)
            << beta * (glob_av[ix]/(double)iblk - pow(glob_av[im]/(double)iblk, 2))/(double)nspin << setw(wd)
            <<  beta * sqrt ( pow(err_x, 2) + pow(2 * fabs(glob_av[im])/(double)iblk * err_m, 2) )/(double)nspin << endl ;

    Chi.close() ;
   
    if (iblk==nblk)
    {
        ofstream Final ;
        if (metro == 1)
            Final.open("Metro/output.final_M" + to_string(h), ios::app) ;
        else
            Final.open ("Gibbs/output.final_M" + to_string(h), ios::app) ;
        Final << temp << setw(wd)
              << glob_av[iu]/(double)iblk/(double)nspin << setw(wd) //ene
              << err_u/(double)nspin << setw(wd)  //ene_err
              << beta*beta * (glob_av[ic]/(double)iblk - pow (glob_av[iu]/(double)iblk, 2))/(double)nspin << setw(wd) //heat
              << beta*beta * sqrt ( pow ( err_c, 2 ) + pow ( 2 * fabs(glob_av[iu])/(double)iblk * err_u, 2 ) )/(double)nspin << setw(wd) //heat_err
              << glob_av[im]/(double)iblk/(double)nspin << setw(wd) //mag
              << err_m/(double)nspin << setw(wd) //mag_err
              << beta * glob_av[ix]/(double)iblk/(double)nspin  << setw(wd) //chi
              << beta * err_x /(double)nspin //chi_err
        << endl;
        
    }
    
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
    ofstream WriteConf;

    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");
    for (int i=0; i<nspin; ++i)
    {
        WriteConf << s[i] << endl;
    }
    WriteConf.close();

    rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
