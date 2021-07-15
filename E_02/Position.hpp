//  Class to describe a position in a 3D world

#ifndef Position_hpp
#define Position_hpp

#include <iostream>
#include <cmath>



class Position
{
private:
    double m_x, m_y, m_z ;
    
protected:
    
public:
    Position () ; //Constructor with no arguments: initialises the data members as 0
    Position ( double x, double y, double z ) ; //Constructor starting from the three cartesian coordinates
    ~Position () {} ;
    
    double x() const ; //Returns m_x
    double y()  const;  //Returns m_y
    double z() const ;   //Returns m_z
    
    void Set_x ( double ) ;  //Sets m_x
    void Set_y ( double ) ;   //Sets m_y
    void Set_z ( double ) ;   //Sets m_z
    
    double r () const ;  //Returns the distance from the origin
    double r2 () const ;  //Returns the square of the distance from the origin
    
    void StepCart ( double stepx, double stepy, double stepz ) ; //Changes the data members according to the step given (in Cartesian coordinates)
    void StepSph ( double phi, double teta, double step) ;  //Changes the data members according to the step given (in spherical coordinates)
    
} ;

#endif /* Position_hpp */
