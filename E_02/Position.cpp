#include "Position.hpp"

Position::Position ()
{
    m_x = 0. ;
    m_y = 0. ;
    m_z = 0. ;
}

Position::Position ( double x, double y, double z )
{
    m_x = x ;
    m_y = y ;
    m_z = z ;
}


double Position::x() const
{
    return m_x ;
}
    
double Position::y()  const
{
    return m_y ;
}

double Position::z() const
{
    return m_z ;
}
   
void Position::Set_x ( double val )
{
    m_x = val;
}
void Position::Set_y ( double val )
{
    m_y = val;
}
void Position::Set_z ( double val )
{
    m_z = val;
}
   
double Position::r () const
{
    return sqrt ( (m_x * m_x) + (m_y * m_y) + (m_z * m_z) ) ;
}

void Position::StepCart ( double stepx, double stepy, double stepz )
{
    m_x += stepx ;
    m_y += stepy ;
    m_z += stepz ;
}

 double Position::r2 () const
 {
     return (m_x * m_x) + (m_y * m_y) + (m_z * m_z) ;
 }

void Position::StepSph ( double phi, double teta, double step)
{
    m_x += sin (teta) * cos (phi) ;
    m_y += sin (teta) * sin (phi) ;
    m_z += cos (teta) ;
}
