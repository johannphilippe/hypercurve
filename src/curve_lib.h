#ifndef CURVE_LIB_H
#define CURVE_LIB_H

#include"utilities.h"
#include<cmath>
#include<iostream>
#include<functional>
#include<vector>
#include"cubic_spline.h"

namespace  hypercurve {

//////////////////////////////////////////////////
// Simple Curves
//////////////////////////////////////////////////
class curve_base
{
public:
    curve_base() {}
    virtual ~curve_base() {}

    // This method is the one to implement/override for simple curves (who only need x and constants to calculate y)
    inline virtual double process(double x) {return x;};
    // When you want to retrieve a single sample from a curve, it is recommanded to use this one. It is not necessary to override it for simple curves,
    // but may be necessary for complex (bezier, spline, catmullrom ...)
    inline virtual double process(size_t i, size_t size) {return process(frac(i, size));}
    inline virtual double process_all(size_t size, std::vector<double>::iterator &it)
    {
        double max = 0.0;
        for(size_t i = 0; i < size; ++i)
        {
            double res = scale(process(frac(i, size)));
            if( std::abs(res) > max) max = std::abs(res);
            *it = res;
            ++it;
        }
        return max;
    }

    inline virtual void init(double y_start_, double y_dest_) {
        y_start = y_start_;
        y_destination = y_dest_;
        abs_diff = std::abs(y_start - y_destination);
        offset = std::min(y_start, y_destination);
        on_init();
    };

    // Override this one insted of init to avoid y_start and y_destination affectation repetition
    inline virtual void on_init() {}

protected:

    double scale(double y)
    {
        if(y_start > y_destination) y = 1.0 - y;
        return (y * abs_diff) + offset;
    }

    double y_start, y_destination, abs_diff, offset;
};

using linear_curve = curve_base;

class diocles_curve : public curve_base
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

class cubic_curve: public curve_base
{
public:
    inline double process(double x) override
    {
        return std::pow(x, 3);
    }
};

// Choose the exponent of X
class exponent_curve  : public curve_base
{
public:
    exponent_curve(double exponent_)
        : exponent(exponent_)
    {}
    inline double process(double x) override
    {
        double res = std::pow(x, exponent) ;
        return res;
    }

private:
    double exponent;
};

//////////////////////////////////////////////////
// Chebyshev
// -- pretty experimental, no idea to make it go from one point(y) to another(y)
// -- could be a modulator instead
//////////////////////////////////////////////////

template<int T = 1>
class chebyshev_curve : public curve_base
{
public:
    chebyshev_curve(int n_)
        : n(n_)
    {}


    inline double process(double x) override
    {

        const double t = std::acos(( x * 2.0) - 1.0 );
        if constexpr(T == 1)
        {
            return scale_chebyshev(std::cos(n * t));
        } else // T == 2
        {
            return scale_chebyshev(std::sin( (n + 1) * t) / std::sin(t));
        }
    }
private:

    // T = 1 stays between -1 and 1
    // T = 2
    double scale_chebyshev(double y)
    {
        return (y) ;
    }

    const double n;
};

//////////////////////////////////////////////////
// Cardioid
// Idea : user could just give a segment in degrees with a rotation offset, and give y_start and y_dest.
//////////////////////////////////////////////////

//////////////////////////////////////////////////
// User defined curve Curves
//
// Callback should return an y value between 0 and 1
// for each x between 0 and 1
//////////////////////////////////////////////////
class user_defined_curve : public curve_base
{
public:
    user_defined_curve() {}
    user_defined_curve(std::function<double(double)> f)
        : callback(f)
    {}

    inline double process(double x) override
    {
        return callback(x);
    }
protected:
    std::function<double(double)> callback;
};


//////////////////////////////////////////////////
// Bezier Curves
//////////////////////////////////////////////////

class bezier_curve_base : public curve_base
{
public:
    inline double process_all(size_t size, std::vector<double>::iterator &it) override
    {
        double max = 0.0;
        size_t cnt = 0;
        std::pair<double, double> r1, r2;
        r1 = process_bezier(double(cnt) / double(size));
        r2 = process_bezier(double(cnt + 1) / double(size));

        for(size_t i = 1; i <= size; i++)
        {
            const double x = double(i) / double(size);
            while(cnt < size)
            {
                if( (r1.first <= x ) && (r2.first >= x))
                    break;
                r1 = process_bezier( double(cnt) / double(size) );
                r2 = process_bezier( double(cnt + 1) / double(size) );
                ++cnt;
            }

            double relative_x = relative_position(r1.first, r2.first, x);
            double linear_interp = linear_interpolation(r1.second, r2.second, relative_x);
            if(std::abs(linear_interp) > max) max = std::abs(linear_interp);
            *it = linear_interp;
            ++it;
        }

        return max;
    }

    inline virtual double process(double x) override
    {
        throw(std::runtime_error("Unimplemented for Bezier curve"));
    }

    // This one should be implemented instead of the above one (if wanted to process single bezier point)
    inline virtual double process(size_t i, size_t size) override
    {
        return process(frac(i, size));
    }
protected:
    inline virtual std::pair<double, double> process_bezier(double x) {return {0,0};}
};

class quadratic_bezier_curve : public bezier_curve_base
{
public:
    quadratic_bezier_curve(control_point cp)
        : _control_point(cp)
    {}

    void on_init() override
    {
        c_x = (3 * (_control_point.x));
        b_x = (-c_x);
        a_x = (1 - c_x - b_x);
        c_y = ( 3 * (_control_point.y - y_start));
        b_y = (-c_y);
        a_y = (y_destination - y_start - c_y - b_y);
    }

private:
    inline std::pair<double, double> process_bezier(double x) override
    {
        return {
            (a_x * std::pow(x, 3)) + (b_x * std::pow(x, 2)) + (c_x * x),
            (a_y * std::pow(x,3)) + (b_y * std::pow(x,2)) + (c_y * x) + y_start
        };
    }

    control_point _control_point;
    double c_x, b_x, a_x, c_y, b_y, a_y;
};

class cubic_bezier_curve : public bezier_curve_base
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

private:
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
    control_point _control_point1, _control_point2;
    double c_x, b_x, a_x, c_y, b_y, a_y;
};

//////////////////////////////////////////////////
// Spline Curves
//////////////////////////////////////////////////

class cubic_spline_curve : public virtual curve_base
{
public:
    cubic_spline_curve(std::vector<point> cp)
        : _control_points( std::move(cp) )
    {
        if(_control_points.size() < 3)
            throw(std::runtime_error("Control point list size must be >= 3"));
    }

    void on_init() override
    {
       _control_points.insert(_control_points.begin(), curve_point(0, y_start));
       //_control_points.push_back(curve_point(1, y_destination));
        // For each point, determine absolute position from relative position
        for(size_t i = 1; i < _control_points.size(); i++)
        {
            _control_points[i].x += _control_points[i-1].x;
        }
    }

    inline virtual double process_all(size_t size, std::vector<double>::iterator &it) override
    {
        std::vector<double>& res = spl.interpolate_from_points(_control_points, size, point{1.0, 1.0});
        double max = 0.0;
        for(size_t i = 0; i < size; i++)
        {
            if( std::abs(res[i]) > max ) max = std::abs(res[i]);
            *it = res[i];
            ++it;
        }
        return max;
    }

private:
    cubic_spline<double> spl;
    std::vector<point> _control_points;
};

// User passes control points P0 and P3 assuming that P1(0,y) and P2(1,y). Calculation will be relative, and will be  rescaled after.
// Based on https://www.desmos.com/calculator/552cpvzfxw?lang=fr

// To get a real centripetal or chordal, we should be based on https://www.desmos.com/calculator/9kazaxavsf?lang=fr
// alpha =  Parametric constant: 0.5 for the centripetal spline, 0.0 for the uniform spline, 1.0 for the chordal spline.
class catmull_rom_spline_curve : public curve_base
{
public:
    catmull_rom_spline_curve(double alpha_, point p0, point p3)
        : alpha(alpha_)
        , _cp0(p0)
        , _cp3(p3)
    {}


    void on_init() override
    {
        _cp1 = point(0, y_start);
        _cp2 = point(1, y_destination);
    }


    inline virtual double process_all(size_t size, std::vector<double>::iterator &it) override
    {
        size_t cnt = 0;
        std::pair<double, double> r1, r2;
        r1 = process_catmul_rom(0);
        r2 = process_catmul_rom(1.0 / double(size));
        double max = 0.0;
        for(size_t i =0; i < size; i++)
        {
            const double x = double(i) / double(size);
            while(cnt < size)
            {
                if((r1.first <= x) && (r2.first >= x))
                    break;
                r1 = process_catmul_rom( double(cnt) / double(size) );
                r2 = process_catmul_rom( double(cnt + 1) / double(size) );
                ++cnt;
            }

            double relative_x = relative_position(r1.first, r2.first, x);
            double linear_interp = linear_interpolation(r1.second, r2.second, relative_x);
            if(std::abs(linear_interp) > max) max = std::abs(linear_interp);
            *it = linear_interp;
            ++it;
        }
        return  max;
    }
private:

    std::pair<double, double> process_catmul_rom(double x)
    {
        const double y = linear_interpolation(y_start, y_destination, x);
        const double rx = 0.5 * ((_cp1.x * 2.0) + (-_cp0.x + _cp2.x) * x
                                 + ((_cp0.x * 2.0) - (_cp1.x * 5.0)
                                    + (_cp2.x * 4.0) - _cp3.x ) * (x*x)
                                 + (-_cp0.x  + (_cp1.x * 3.0) - (_cp2.x * 3.0) +_cp3.x) * (x*x*x));
        const double ry = 0.5 * ((_cp1.y * 2.0) + (-_cp0.y + _cp2.y) * x
                                 + ((_cp0.y * 2.0) - (_cp1.y * 5.0)
                                    + (_cp2.y * 4.0) - _cp3.y ) * (x*x)
                                 + (-_cp0.y  + (_cp1.y * 3.0) - (_cp2.y * 3.0) +_cp3.y) * (x*x*x));
        return {rx, ry};
    }

    double alpha = 0;
    control_point _cp0, _cp3, _cp1, _cp2;
    double t0, t1, t2, t3;
};


}
#endif // CURVE_LIB_H

