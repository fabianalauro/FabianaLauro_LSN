/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

class Random
{
private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random();     // Constructor with no arguments
  Random (std::string, std::string) ;     //Constructor starting from the name of the file where p1 and p2 are written and the name of the file where the 4 seeds are written, preceded by "RANDOMSEED"
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);    //Sets the generator with an array with the 4 seeds and two integers p1 e p2
  void SetRandom (std::string, std::string) ;   //Sets the generator with the name of the file where p1 and p2 are written and the name of the file where the 4 seeds are written, preceded by "RANDOMSEED"
  void SaveSeed();   //Prints on the file "seed.out" l1, l2, l3 e l4
  double Rannyu(void);       //Generates random numbers uniformly distributed between 0 e 1
  double Rannyu(double min, double max);      //Generates random numbers uniformly distributed between min and max
  double Gauss(double mean, double sigma);      //Generates random numbers with Gaussian distribution
  double Exp (double lambda) ;      //Generates random numbers with exponential distribution
  double Lorentz (double x0, double gamma) ; //Generates random numbers with Cauchy-Lorentz distribution
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
