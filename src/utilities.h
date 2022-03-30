#ifndef UTILITIES_H
#define UTILITIES_H

#include<math.h>
#include<iostream>

namespace hypercurve {

struct point
{
    point(double x_, double y_)
        : x(x_)
        , y(y_)
    {}
    double x,y;
};

// x must be between 0 and 1, 0 being the y1 x position, 1 the y2 x position
inline double linear_interpolation(double y1, double y2, double x)
{
    return y1 + (x * (y2 - y1)) ;
}

// returns x between 0 and 1 (1 being x == x2, 0 being x == x1)
// make sure x1 <= x <= x2
inline double relative_position(double x1, double x2, double x)
{
    if(!(x >= x1 && x <= x2) || !(x1 < x2))
        throw(std::runtime_error("Make sure x1 <= x <= x2 and x1 < x2"));
    const double factor = 1.0 / (x2 - x1);
    return (x * factor) - (x1 * factor);
}

}

#endif // UTILITIES_H
