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
using control_point = point;

struct curve_point : public point
{
    curve_point(float x_, float y_, float curve_ = 0.0f) :
        point(x_, y_), curve(curve_)
    {}

    curve_point(point &p, float curve_ = 0.0f) :
        point(p), curve(curve_)
    {}

    template<typename T>
    float distance_to(T &p)
    {
        const float a2 = (p.x - x) * (p.x - x);
        const float b2 = (p.y - y) * (p.y - y);
        const float dist = sqrt(a2 + b2);
        return dist;
    }

   float curve;
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

inline double limit(double min, double max, double v)
{
    if(v > max) return max;
    if(v < min) return min;
    return v;
}

}

#endif // UTILITIES_H
