
#include "Probability.hpp"


double Hydr100::Probability (double x, double y, double z) 
{
    return 1./M_PI * exp ( -2. * sqrt (x*x + y*y + z*z) ) ;
}

double Hydr210::Probability (double x, double y, double z)
{
    return 1./(32 * M_PI) * z*z * exp ( - sqrt (x*x + y*y + z*z) ) ;
}
