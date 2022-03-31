#ifndef UTILITIES_H
#define UTILITIES_H

#include<math.h>
#include<iostream>
#include<vector>
#include<memory>

namespace hypercurve {


double pos(double x)
{
    return (x > 0) ? x : 0;
}

double frac(double x1, double x2)
{
    return double(x1) / double(x2);
}

template<typename T>
std::shared_ptr<T> share(T t)
{
    return std::make_shared<T>(t);
}

struct point
{
    point(double x_, double y_)
        : x(x_)
        , y(y_)
    {}

    point() {}

    void print()
    {
        printf("point \n\tx = %f \n\ty = %f\n", x, y);
    }

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

// Transforms the results of an interpolation to x relative values (useful for Bezier and Catmul Rom for example)
template<typename T>
static double remap_coordinates(std::vector< std::pair<T, T> > &map, std::vector<T> &res)
{
    size_t cnt = 0;
    double max = -1;
    for(size_t i = 0; i < map.size() - 1; i++)
    {
        const double x = double(i) / double(map.size());
        while(cnt < map.size())
        {
            if((map[cnt].first <= x) && (map[cnt + 1].first >= x))
                break;
            ++cnt;
        }
        const double relative_x = relative_position(map[cnt].first, map[cnt+1].first, x);
        const double linear_interp = linear_interpolation(map[cnt].second, map[cnt+1].second, relative_x);
        if(linear_interp > max) max = linear_interp;
        res[i] = linear_interp;
    }
    return max;
}

std::vector<double> linspace(size_t size)
{
    std::vector<double> v(size);
    for(size_t i = 0; i < size; i++)
        v[i] = (double(i) / double(size));
    return v;
}

}

#endif // UTILITIES_H
