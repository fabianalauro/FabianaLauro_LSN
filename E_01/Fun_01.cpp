
#include <iostream>
#include <cmath>
#include "Fun_01.h"

double dist (double x1, double y1, double x2, double y2)
{
    return std::sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
}

double slope (double x1, double y1, double x2, double y2)
{
    return (y2-y1) / (x2-x1) ;
}

double y_cf ( double m, double r, bool positive )
{
    if (positive)
        return std::abs(m) * r / std::sqrt( m*m + 1. ) ;
    else
        return - std::abs(m) * r / std::sqrt( m*m + 1. ) ;
}
