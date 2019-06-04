#ifndef _POINT_HPP
#define _POINT_HPP

#include "../include/point.h"

Point::Point()
{
    // Default Point constructor
    // Sets x = 0 and y = 0
    x = 0;
    y = 0;
}

Point::Point(const double xx, const double yy)
{
    // Constructs Point object (x = xx, y = yy)
    x = xx;
    y = yy;
}

void Point::setX(const double xx)
{
    // Sets x = xx
    x = xx;
}

void Point::setY(const double yy)
{
    // Sets y = yy
    y = yy;
}

double Point::getX()
{
    // Gets x = xx
    return x;
}

double Point::getY()
{
    // Gets y = yy
    return y;
}
;
#endif