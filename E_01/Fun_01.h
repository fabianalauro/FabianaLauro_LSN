//Functions created for Exercise 1

#ifndef Fun_01_h
#define Fun_01_h

double dist ( double x1, double y1, double x2, double y2 ) ; //Returns distance between two points

double slope ( double x1, double y1, double x2, double y2 ) ;  //returns the slope of the line between two points

double y_cf ( double m, double r, bool positive ) ;  //Returns the y-coordinate of the intersection between a line and a circunference, both centered in the origin of the axes. The bool is needed to choose between the two solutions

#endif /* Fun_01_hpp */
