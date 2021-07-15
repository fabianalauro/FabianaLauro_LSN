/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Molecular dynamics simulation in a NVE ensemble, with a classic Lennard-Jones fluid, using the Verlet method

#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow, sqrt
#include <string>       // std::string
#include "MolDyn_NVE.h"

using namespace std;

int main (int argc, char ** argv)
{
    if ( argc != 4 )
    {
        cout << "Usage: <" << argv [0] << "> <readold> <rescale> <average>" << endl << endl ;
        exit (1) ;
    }
   
    readold = atoi( argv[1] ) ; //Whether to read from file old.0
    rescale = atoi ( argv[2] ) ;  //Whether to rescale velocities
    average = atoi (argv[3]) ;  //Whether to compute average values
        
    Input();             //Initialisation
    int nconf = 1;
    
    //Initialising classes for data blocking
    BA_ekin.Set_parameters (nstep/10, nblocks) ;
    BA_epot.Set_parameters (nstep/10, nblocks) ;
    BA_etot.Set_parameters (nstep/10, nblocks) ;
    BA_temp.Set_parameters (nstep/10, nblocks) ;
    for (int i = 0; i < m_props; i++)
        BA_gofr[i].Set_parameters (nstep/10, nblocks) ;
    
    //Simulation
    for(int istep = 1; istep <= nstep; ++istep)
    {
        for (int i = 0; i < nbins; i++)
            bins[i] = 0 ;
        Move();           //Move particles with Verlet algorithm

        if(istep%iprint == 0)
        {
            cout << "Number of time-steps: " << istep << endl;
        }
        if(istep%10 == 0)
        {
            Measure();     //Properties measurement
            ConfXYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
            nconf += 1;
        }
    }
    ConfFinal();         //Write final configuration to restart

    return 0;
}


//Prepare all stuff for the simulation
void Input(void)
{
    ifstream ReadInput,ReadConf;
    double ep, ek, pr, et, vir;

    cout << "Classic Lennard-Jones fluid        " << endl;
    cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
    cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
    cout << "The program uses Lennard-Jones units " << endl;

    seed = 1;    //Set seed for random numbers
    srand(seed); //Initialize random number generator
  
    ReadInput.open("input.dat"); //Read input

    ReadInput >> temp;

    ReadInput >> npart;
    cout << "Number of particles = " << npart << endl;

    ReadInput >> rho;
    cout << "Density of particles = " << rho << endl;
    vol = (double)npart/rho;
    cout << "Volume of the simulation box = " << vol << endl;
    box = pow(vol,1.0/3.0);
    cout << "Edge of the simulation box = " << box << endl;

    ReadInput >> rcut;
    ReadInput >> delta;
    ReadInput >> nstep;
    ReadInput >> nblocks;
    ReadInput >> iprint;

    cout << "The program integrates Newton equations with the Verlet method " << endl;
    cout << "Time step = " << delta << endl;
    cout << "Number of steps = " << nstep << endl << endl;
    if (average == 1)
    {
        cout << "Computing averages using data blocking" << endl ;
        cout << "Number of blocks = " << nblocks << endl << endl;
    }

    ReadInput.close();

    //Prepare array for measurements
    iv = 0; //Potential energy
    ik = 1; //Kinetic energy
    ie = 2; //Total energy
    it = 3; //Temperature
    n_props = 4; //Number of observables
    nbins = 100;
    bin_size = (box/2.0)/(double)nbins;
    
    //Read initial configuration
    cout << "Read initial configuration from file config.0" << endl << endl;
    ReadConf.open("config.0");
    for (int i = 0; i < npart; ++i)
    {
        ReadConf >> x[i] >> y[i] >> z[i];
        x[i] = x[i] * box;
        y[i] = y[i] * box;
        z[i] = z[i] * box;
    }
    ReadConf.close();

//Generating old initial configuration
    if (readold == 0)
    {
       //Prepare initial velocities
       cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
       double sumv[3] = {0.0, 0.0, 0.0};
       for (int i = 0; i < npart; ++i)
       {
           vx[i] = rand()/double(RAND_MAX) - 0.5;
           vy[i] = rand()/double(RAND_MAX) - 0.5;
           vz[i] = rand()/double(RAND_MAX) - 0.5;

           sumv[0] += vx[i];
           sumv[1] += vy[i];
           sumv[2] += vz[i];
       }
       for (int idim = 0; idim < 3; ++idim)
           sumv[idim] /= (double)npart;
       double sumv2 = 0.0 ;
       for (int i = 0; i < npart; ++i)
       {
           vx[i] = vx[i] - sumv[0];
           vy[i] = vy[i] - sumv[1];
           vz[i] = vz[i] - sumv[2];

           sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
       }
       sumv2 /= (double)npart;

        Rescale (sumv2) ;
    }
    
//Import old initial configuration from file and (optional) rescale velocities according to temperature
    
    else if (readold == 1)
    {
        cout << "Importing old initial configuration from file old.0 " << endl  ;
        ReadConf.open("old.0");
        for (int i = 0; i < npart; ++i)
        {
            ReadConf >> xold[i] >> yold[i] >> zold[i];
            xold[i] = xold[i] * box;
            yold[i] = yold[i] * box;
            zold[i] = zold[i] * box;
        }
        ReadConf.close();
        
        if (rescale == 1 )
        {
            cout << "Rescaling velocities according to the given temperature" << endl << endl ;
        
            Move () ;
        
            double sumv2 = 0. ;
        
            for (int i = 0; i < npart; i++)
            {
                vx[i] = (x[i] - xold[i])  / delta ;
                vy[i] = (y[i] - yold[i])  / delta ;
                vz[i] = (z[i] - zold[i])  / delta ;
                sumv2 += pow( (x[i] - xold[i])  / delta, 2 ) + pow( (y[i] - yold[i]) / delta, 2 ) + pow( (z[i] - zold[i]) / delta, 2 ) ;
            }
            sumv2 /= (double)npart ;
            Rescale (sumv2) ;
        }
    }
    
   return;
}

void Rescale (double vqm)  //Initialise old coordinates rescaling the velocities according to the temperature
{
    double fs = sqrt(3 * temp / vqm);   // fs = velocity scale factor (temp/T)
    //cout << fs ;
    for (int i = 0; i < npart; ++i)
    {
        vx[i] *= fs;
        vy[i] *= fs;
        vz[i] *= fs;
        
        xold[i] = Pbc(x[i] - vx[i] * delta);
        yold[i] = Pbc(y[i] - vy[i] * delta);
        zold[i] = Pbc(z[i] - vz[i] * delta);
    }
    return ;
}

void Move(void) //Move particles with Verlet algorithm
{
    double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

    for(int i = 0; i < npart; ++i) //Force acting on particle i
    {
        fx[i] = Force(i,0);
        fy[i] = Force(i,1);
        fz[i] = Force(i,2);
    }

    for(int i = 0; i < npart; ++i) //Verlet integration scheme
    {
        xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
        ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
        znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
        
        vx[i] = Pbc(xnew - xold[i]) / (2.0 * delta);
        vy[i] = Pbc(ynew - yold[i]) / (2.0 * delta);
        vz[i] = Pbc(znew - zold[i]) / (2.0 * delta);
        
        xold[i] = x[i];
        yold[i] = y[i];
        zold[i] = z[i];
        
        x[i] = xnew;
        y[i] = ynew;
        z[i] = znew;
    }
    return;
}

double Force(int ip, int idir) //Compute forces as -Grad_ip V(r)
{
    double f = 0.0;
    double dvec[3], dr;
    
    for (int i = 0; i < npart; ++i)
    {
        if(i != ip)
        {
            dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
            dvec[1] = Pbc( y[ip] - y[i] );
            dvec[2] = Pbc( z[ip] - z[i] );
            
            dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
            dr = sqrt(dr);
            
            if(dr < rcut)
            {
                f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
            }
        }
    }
    
    return f;
}

void Measure() //Properties measurement
{
    int bin ;
    double v, t, vij;
    double dx, dy, dz, dr;
    ofstream Epot, Ekin, Etot, Temp, gave;
    
    v = 0.0; //reset observables
    t = 0.0;
    
    //cycle over pairs of particles
    for (int i = 0; i < npart-1; ++i)
    {
        for (int j = i+1; j < npart; ++j)
        {
            
            dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
            dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
            dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)
            
            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
            
            bin = static_cast<int>( trunc(dr / bin_size) );
            bins[bin] += 2 ;
            
            
            
            /*************************
            cout << "DISTANCE : g(r)" << endl;
               for (int k=0; k<nbins; ++k)
               {
                   double r = k * bin_size;
                   double DeltaV = 4.*M_PI/3. * ( pow(r+bin_size, 3) - pow (r, 3) ) ;
                   cout << bins [k] / (rho*npart*DeltaV) << " " ;
               }
               cout << endl;
             ************************/
           
            
            
            if(dr < rcut)
            {
                vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
                
                //Potential energy
                v += vij;
            }
        }
    }
    
    //Kinetic energy
    for (int i = 0; i < npart; ++i)
        t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
    //cout <<t ;
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    //Opening the right files and writing on them
    if (average == 0)
    {
        if ( temp == 0.8 )
        {
            Epot.open("S_output_epot.dat",ios::app);
            Ekin.open("S_output_ekin.dat",ios::app);
            Temp.open("S_output_temp.dat",ios::app);
            Etot.open("S_output_etot.dat",ios::app);
        }
    
        else if ( temp == 1.1 )
        {
            Epot.open("L_output_epot.dat",ios::app);
            Ekin.open("L_output_ekin.dat",ios::app);
            Temp.open("L_output_temp.dat",ios::app);
            Etot.open("L_output_etot.dat",ios::app);
        }
      
        else if ( temp == 1.2 )
        {
            Epot.open("G_output_epot.dat",ios::app);
            Ekin.open("G_output_ekin.dat",ios::app);
            Temp.open("G_output_temp.dat",ios::app);
            Etot.open("G_output_etot.dat",ios::app);
        }
        Epot << stima_pot  << endl;
        Ekin << stima_kin  << endl;
        Temp << stima_temp << endl;
        Etot << stima_etot << endl;

        Epot.close();
        Ekin.close();
        Temp.close();
        Etot.close();
    }
    
    if (average == 1)
    {
        if ( temp == 0.8 )
        {
            Epot.open ("S_ave_epot.out", ios::app);
            Ekin.open ("S_ave_ekin.out", ios::app);
            Etot.open ("S_ave_etot.out", ios::app);
            Temp.open ("S_ave_temp.out", ios::app);
            gave.open ("S_ave_gofr.dat", ios::app);
        }
        else if ( temp == 1.1 )
        {
            Epot.open ("L_ave_epot.out", ios::app);
            Ekin.open ("L_ave_ekin.out", ios::app);
            Etot.open ("L_ave_etot.out", ios::app);
            Temp.open ("L_ave_temp.out", ios::app);
            gave.open ("L_ave_gofr.dat", ios::app);
        }
        else if ( temp == 1.2 )
        {
            Epot.open ("G_ave_epot.out", ios::app);
            Ekin.open ("G_ave_ekin.out", ios::app);
            Etot.open ("G_ave_etot.out", ios::app);
            Temp.open ("G_ave_temp.out", ios::app);
            gave.open ("G_ave_gofr.dat", ios::app);
        }
        
        BA_epot.Add(stima_pot, Epot) ;
        BA_ekin.Add(stima_kin, Ekin) ;
        BA_temp.Add(stima_temp, Temp) ;
        BA_etot.Add(stima_etot, Etot) ;
        
        for (int i=0; i<nbins; i++)
        {
            double r = i * bin_size;
            double DeltaV = 4.*M_PI/3. * ( pow(r+bin_size, 3) - pow (r, 3) ) ;
            double stima_gdir = bins[i] / (rho*npart*DeltaV) ;
            BA_gofr[i].Writefinal ( stima_gdir, nblocks, gave ) ;
        }
        
        Epot.close();
        Ekin.close();
        Temp.close();
        Etot.close();
        gave.close() ;
        
    }
    
    return;
}


void ConfFinal(void) //Write final configuration and the final old configuration in two different files
{
    ofstream WriteConf;

    cout << "Print final configuration to file config.final " << endl << endl;
    WriteConf.open("config.final");

    for (int i = 0; i < npart; ++i)
    {
        WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
    }
    WriteConf.close();
    
    cout << "Print final old configuration to file old.final " << endl << endl;
    WriteConf.open ("old.final") ;
    
    for (int i = 0; i < npart; ++i)
       {
           WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
       }
    WriteConf.close();
    
    return;
}

void ConfXYZ(int nconf) //Write configuration in .xyz format
{
    ofstream WriteXYZ;

    WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
    WriteXYZ << npart << endl;
    WriteXYZ << "This is only a comment!" << endl;
    for (int i = 0; i < npart; ++i)
    {
        WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
    }
    WriteXYZ.close();
}

double Pbc(double r) //Algorithm for periodic boundary conditions with side L=box
{
    return r - box * rint(r/box);
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
