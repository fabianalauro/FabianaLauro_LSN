//Classes to implement various probability distribution functions, using polimorphism
//CLASSES USED BY THE CLASS "METROPOLIS"

#ifndef Probability_hpp
#define Probability_hpp

#include <cmath>

class Basic_PDF
{
private:
    
protected:
    

public:
    
    Basic_PDF () {} ;
    virtual ~Basic_PDF () {} ;
    
    virtual double Probability (double x, double y, double z) = 0 ;
    
} ;


class Hydr100 : public Basic_PDF
{
    private:
        
    protected:
        
    public:
    
    Hydr100 () {} ;
    ~Hydr100 () {} ;
    
    double Probability (double x, double y, double z) ;
    
};


class Hydr210 : public Basic_PDF
{
    private:
        
    protected:
        
    public:
    
     ~Hydr210 () {} ;
    
    double Probability (double x, double y, double z) ;

};

class Variational : public Basic_PDF
{
private:
    double m_sigma, m_mu ;
protected:
    
public:
    
    Variational () ;
    Variational (double sigma, double mu) ;
    ~Variational () {} ;
    
    void set_parameters (double sigma, double mu) ;
    
    double Probability (double x, double y, double z) ;
} ;



#endif /* Probability_hpp */
