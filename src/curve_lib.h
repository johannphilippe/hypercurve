#ifndef CURVE_LIB_H
#define CURVE_LIB_H

#include"utilities.h"
#include<math.h>
#include<iostream>
#include<functional>

namespace  hypercurve {

//////////////////////////////////////////////////
// Simple Curves
//////////////////////////////////////////////////
class curve_base
{
public:
    curve_base() {}
    virtual ~curve_base() {}

    inline virtual double process(double x) {return x;};

    inline virtual void init(double y_start, double y_dest) {};
};

using linear_curve = curve_base;

class diocles_curve : public virtual curve_base
{
public:
    diocles_curve(const double a_)
        : a(a_)
    {}
    inline double process(double x) override
    {
        return std::sqrt( std::pow(x, 3) / (2 * a - x) );
    }

private:
    const double a;
};

using cissoid_curve = diocles_curve;

class cubic_curve: public virtual curve_base
{
public:
    inline double process(double x) override
    {
        return std::pow(x, 3);
    }
};

// Custom homemade curve
class hypersmooth_curve  : public curve_base
{
public:
    inline double process(double x) override
    {
        return std::pow(x, 3) - std::pow(x / 2.0, 2) * log(std::pow(x, 4));
    }
};

//////////////////////////////////////////////////
// User defined curve Curves
//////////////////////////////////////////////////
class user_defined_curve : public curve_base
{
public:
    user_defined_curve(std::function<double(double)> f)
        : callback(f)
    {}

    inline double process(double x) override
    {
        return callback(x);
    }
private:
    std::function<double(double)> callback;
};


//////////////////////////////////////////////////
// Bezier Curves
//////////////////////////////////////////////////

using control_point = point;

class bezier_curve_base : public virtual curve_base
{
public:
    void init(double y_start, double y_dest) override
    {}

    virtual std::pair<double, double> process_bezier(double x) {
        return {0,0};
    }

protected:
    double y_start, y_destination;
};

class quadratic_bezier_curve : public virtual bezier_curve_base
{
public:
    quadratic_bezier_curve(control_point cp)
        : _control_point(cp)
    {}

    void init(double y_start_, double y_dest_) override
    {
        y_start = y_start_;
        y_destination = y_dest_;
        c_x = (3 * (_control_point.x));
        b_x = (-c_x);
        a_x = (1 - c_x - b_x);
        c_y = ( 3 * (_control_point.y - y_start_));
        b_y = (-c_y);
        a_y = (y_dest_ - y_start_ - c_y - b_y);
    }

    // Not used in the library, but this might be closer to canonical form since it calculates x as well
    inline std::pair<double, double> process_bezier(double x) override
    {
        return {
            (a_x * std::pow(x, 3)) + (b_x * std::pow(x, 2)) + (c_x * x),
            (a_y * std::pow(x,3)) + (b_y * std::pow(x,2)) + (c_y * x) + y_start
        };
    }

    inline double process(double x) override
    {
        return (a_y * std::pow(x,3)) + (b_y * std::pow(x,2)) + (c_y * x) + y_start;
    }

private:
    control_point _control_point;
    double c_x, b_x, a_x, c_y, b_y, a_y;
};

class cubic_bezier_curve : public virtual bezier_curve_base
{
public:
    cubic_bezier_curve(control_point _cp1, control_point _cp2)
        : _control_point1(_cp1)
        , _control_point2(_cp2)
    {}

    void init(double y_start_, double y_dest_) override
    {
        y_start = y_start_;
        y_destination = y_dest_;
        c_x = ((-1) * 0 + 3 * _control_point1.x - 3 * _control_point2.x + 1);
        b_x = (3 * 0 - 6 * _control_point1.x + 3 * _control_point2.x);
        a_x = ((-3) * 0 + 3 * _control_point1.x);
        c_y = ((-1) * y_start + 3 * _control_point1.y - 3 * _control_point2.y + y_destination);
        b_y = (3 * y_start - 6 * _control_point1.y + 3 * _control_point2.y);
        a_y = ((-3) * y_start + 3 * _control_point1.y);
    }

    std::pair<double, double> process_bezier(double x) override
    {
        return {
            (std::pow(x, 3) * c_x)
                    + (std::pow(x, 2) * b_x)
                    + (x * a_x),
            (std::pow(x, 3) * c_y)
                    + (std::pow(x, 2) * b_y)
                    + (x * a_y)
                    + y_start
        };
    }
    inline double process(double x) override
    {
            return (std::pow(x, 3) * c_y)
                    + (std::pow(x, 2) * b_y)
                    + (x * a_y)
                    + y_start;
    }

private:
    control_point _control_point1, _control_point2;
    double c_x, b_x, a_x, c_y, b_y, a_y;
};

}
#endif // CURVE_LIB_H
