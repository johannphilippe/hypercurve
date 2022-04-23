/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

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
// Base Curve
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
    inline virtual double process_all(size_t size, memory_vector<double>::iterator &it)
    {
        double max = 0.0;
        for(size_t i = 0; i < size; ++i)
        {
            const double x = frac(i, size);
            double res = scale(x);
            if( std::abs(res) > max) max = std::abs(res);
            if(inverted) res = process_invert(x, res);
            *it = res;
            ++it;
        }
        return max;
    }

    // Do not override this (or make sure to implement members definition)
    inline virtual void init(double y_start_, double y_dest_, size_t definition_) {
        y_start = y_start_;
        y_destination = y_dest_;
        definition = definition_;
        abs_diff = std::abs(y_start - y_destination);
        offset = std::min(y_start, y_destination);
        on_init();
    };

    // Override this one insted of init to avoid y_start and y_destination
    // affectation repetition
    inline virtual void on_init() {}

    // Allows inversion of curve (y symetry on a linear x_start/x_end axis).
    bool inverted = false;
protected:

    inline virtual double scale(double x)
    {
        if(y_start > y_destination)
            return  process(1.0 - x) * abs_diff + offset;
        return process(x) * abs_diff + offset;
    }

    inline double process_invert(double x, double y)
    {
        const double lin = linear_interpolation(y_start, y_destination, x);
        return lin + (lin - y);
    }

    size_t definition;
    double y_start, y_destination, abs_diff, offset;
};

using linear_curve = curve_base;

// Allows you to make an y symetry on a x_start/y_start x_end/y_end linear axis
inline std::shared_ptr<curve_base> invert(std::shared_ptr<curve_base> cb)
{
    cb->inverted = true;
    return cb;
}

//////////////////////////////////////////////////
// Curve Library
//////////////////////////////////////////////////

class diocles_curve : public curve_base
{
public:
    diocles_curve(double a_)
        : a(a_)
    {}
    inline double process(double x) override
    {
        return std::sqrt( std::pow(x, 3) / (2 * a - x) );
    }

private:
    double a;
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
class power_curve  : public curve_base
{
public:
    power_curve(double exponent_)
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

class hanning_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        return  hanning(x * definition, definition * 2);
    }
};

class hamming_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        return  hamming(x * definition, definition * 2);
    }
};


class blackman_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        return  blackman(x * definition, definition * 2);
    }
};

// See https://mathcurve.com/courbes2d/bouche/bouche.shtml
class mouse_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        const double _x = x -1;
        return std::sqrt( std::pow( (a*a) - (_x*_x), 3 ) / std::pow(a, 4) );
    }

    constexpr static const double a = 1.0;
};
using kiss_curve = mouse_curve;

// See https://mathcurve.com/courbes2d/bicorne/bicorne.shtml
class bicorn_curve : public curve_base
{
public:
    // Sign true = positive, false = negative
    bicorn_curve(bool sign_)
        : sign(sign_ == true ? 1 : -1)
    {}

    inline double process(double x) override
    {
        const double _x = x-1;
        if(sign == 1)
            return process_bicorn(_x) / process_bicorn(0);
        return process_bicorn(_x);
    }

    constexpr static const double a = 1.0;
    int sign;

protected:
    inline double process_bicorn(double x)
    {
        return  ((a*a) - (x*x)) / (2*a + (std::sqrt(a*a - x*x) * sign));
    }
};
using cocked_hat_curve = bicorn_curve;

//////////////////////////////////////////////////
// Typed curve - inspired by Csound GEN 16
//////////////////////////////////////////////////
class typed_curve : public curve_base
{
public:
    typed_curve(double type_)
        : type(-type_)
    {}

    inline double process(double x) override
    {
        return log_exp_point<double>(0, 1, definition, x * definition, type);
    }

    double type;
};

//////////////////////////////////////////////////
// Gauss Curve
//////////////////////////////////////////////////
class gauss_curve :  public curve_base
{
public:
    gauss_curve(const double A_, const double c_)
        : A(A_)
        , c(c_)
    {
        process_width();
    }

    inline double process(double x) override
    {
        const double _x = (x * half_width) - half_width;
        double gauss = A * std::exp(-std::pow(_x, 2) / (2 * (c*c)));
        gauss -= y_offset;
        gauss /= (A - y_offset); //  scale between offset and 1
        return gauss;
    }

protected:

    void process_width()
    {
        half_width = std::sqrt(2 * std::log(10) * c);
        y_offset = A * std::exp(-std::pow(half_width, 2) / (2 * (c*c)));
    }

    double A, c;
    double half_width, y_offset;

};
using gaussian_curve = gauss_curve;


//////////////////////////////////////////////////
// Catenary (funicular) curve
// In french "courbe de chainette"
// https://mathcurve.com/courbes2d.gb/chainette/chainette.shtml
// The more 'a' is important, the more this curve is straight
//////////////////////////////////////////////////
class catenary_curve : public curve_base
{
public:
    catenary_curve(double a_)
        : a(a_)
        , height(catenary_process(1))
    {
        if(!(a>0))
            throw (std::runtime_error("Catenary curve : a must be > 0"));
    }
    double process(double x) override
    {
        return catenary_process(x) / height;
    }
protected:
    double catenary_process(double x)
    {
        return a * (std::cosh(x/a)) - a;
    }

    double a, height;

    constexpr static const double mult = 12;
};
using funicular_curve = catenary_curve;


//////////////////////////////////////////////////
// Conchal curve
// https://mathcurve.com/courbes2d.gb/conchale/conchale.shtml
// TODO invert x y axis
// x must be scaled from -a to max_x
// Must be scaled and axis inverted
// Still to do
// Angle b = -4.713
//////////////////////////////////////////////////
class conchal_curve : public curve_base
{
public:
    conchal_curve(double a_, double c_)
        : a(a_)
        , c(c_)
        , w(std::sqrt((a*a) + (c*c)))
        , height(a + w)
    {
        if( (a > c) || (a <= 0) || (c <= 0) || (c <  (a * 1.2)))
            throw(std::runtime_error("Conchal curve : a must be < c, c must be >= a * 1.2, a and c must be > 0"));
        max_x = a + (std::sqrt((a*a) + (c*c)));
    }

    double process(double x) override
    {
        const double step = height / double(this->definition);
        const double mx = x * height - a + step;
        const double my = conchal_process(mx) / conchal_process(-a);
        std::cout << "mx : " << mx << " my : " << my << std::endl;
        const double res = (rotate(mx, my, pi_angle(-90), [&](double x_) {
           return this->conchal_process(x_);
        }) + a ) / (height);

        //std::cout << "res = " << res << std::endl;
        if(std::isnan(res))
            std::cout <<"NAN ALERT at " << x  << "  mx : " << mx << std::endl;
        return res;
    }

protected:
    inline double conchal_process(double x)
    {
        std::cout <<"conchal : " << ((c*c) + (a*a) - (x*x) ) * ((c*c) - (a*a) + (x*x)) << "   " << (x+a) << std::endl;
        return std::sqrt(
                    ((c*c) + (a*a) - (x*x))
                    * ((c*c) - (a*a) + (x*x))
                    ) / (x+a);
    }

    inline double scaled_conchal(double x)
    {
        double _x = (x * (max_x + a)) - a;
        double _conchal_res = conchal_process(_x);
        return x + (_conchal_res + a);
    }

    double a, c, max_x, w, height;
};

//////////////////////////////////////////////////
// Tightrope walker cure
// https://mathcurve.com/courbes2d.gb/danseur/danseur.shtml
// The more a is sperior to b, the more is is "straight"
// abs(b) is the curve x max point
// Also a must not be too close to b (undefined behavior) : a = 1.01 b = 1 is ok. a = 1.0001 is not.
//////////////////////////////////////////////////

class tightrope_walker_curve : public curve_base
{
public:
    tightrope_walker_curve(double a_, double b_)
        : a(a_)
        , b( std::abs(b_) )
    {
        if( (a < 0) || (a <= b) )
            throw(std::runtime_error("a must be > abs(b), a must be > 0"));
    }

    void on_init() override
    {
        max_x = b - b / double(definition);
        max_y = abs(process_tightrope_walker(max_x));
    }

    double process(double x) override
    {
        if( ( x * b) >= max_x) return process_tightrope_walker(max_x) / max_y;
        return process_tightrope_walker(x * b) / max_y;
    }

protected:
    inline double process_tightrope_walker(double x)
    {
        return (x * (a-x)) / std::sqrt((b*b)-(x*x));
    }

    double a,b, max_x, max_y;
};

//////////////////////////////////////////////////
// Toxoid curve
// AKA duplicatrix cubic
// https://mathcurve.com/courbes2d/cubicduplicatrice/cubicduplicatrice.shtml
// a needs to be <= 0 (see constructor)
//////////////////////////////////////////////////


class toxoid_curve : public curve_base
{

public:
    toxoid_curve(double a_)
        : a(-std::abs(a_))
        , height(toxoid_process(1))
    {}

    double process(double x) override
    {
       return toxoid_process(x) / height;
    }

protected:

    inline double toxoid_process(double x)
    {
        return std::sqrt(x * std::pow(x - a/2, 2));
    }

    double a, height;
};
/*
class toxoid_curve : public curve_base
{

public:
    toxoid_curve(double a_)
        : a(a_)
        , height(toxoid_process(a + 1))
    {
        if(a <= 0.0)
           throw(std::runtime_error("a must be > 0"));
    }

    double process(double x) override
    {
        return toxoid_process(x + (a / 2)) / height;
    }

protected:

    inline double toxoid_process(double x)
    {
        //return std::sqrt(x * std::pow(x - a/2, 2));
        return std::sqrt( ( (x*x) * (x - a/2) ) / a  );
    }

    double a, height;
};
*/
using duplicatrix_cubic = toxoid_curve;

//////////////////////////////////////////////////
// user defined curves
//
// callback should return an y value between 0 and 1
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
// vararg polynomial
//
// If you give three parameters a, b, c
// it will return ax^3 + bx^2 + cx
//////////////////////////////////////////////////
class polynomial_curve : public curve_base
{
public:
    polynomial_curve(memory_vector<double> args)
        : constants(args)
    {}

    polynomial_curve(std::vector<double> args)
        : constants(args)
    {}

    polynomial_curve(double *args, size_t size)
        : constants(args, size)
    {}

    double process(double x) override
    {
        double res = 0;
        for(size_t i = 0; i < constants.size(); ++i)
        {
            res += std::pow(x, constants.size() - i) * constants[i];
        }
        // then scale.
        return res;
    }

    memory_vector<double> constants;
};

//////////////////////////////////////////////////
// Bezier Curves
//////////////////////////////////////////////////

class bezier_curve_base : public curve_base
{
public:
    inline double process_all(size_t size, memory_vector<double>::iterator &it) override
    {
        double max = 0.0;
        size_t cnt = 0;
        std::pair<double, double> r1, r2;
        r1 = process_bezier(double(cnt) / double(size));
        r2 = process_bezier(double(cnt + 1) / double(size));

        for(size_t i = 0; i < size ; i++)
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
            if(inverted) linear_interp = process_invert(x, linear_interp);
            *it = linear_interp;
            ++it;
        }

        return max;
    }

    inline virtual double process(double) override
    {
        throw(std::runtime_error("Unimplemented for Bezier curve"));
    }

    // This one should be implemented instead of the above one (if wanted to process single bezier point)
    inline virtual double process(size_t i, size_t size) override
    {
        return process(frac(i, size));
    }
protected:
    inline virtual std::pair<double, double> process_bezier(double x) {return {x,0};}
};


//https://math.stackexchange.com/questions/3280087/can-the-t-of-a-quadratic-bezier-curve-be-found-given-a-x-or-y-using
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
    double get_t(double x)
    {
        return 0;
    }

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


// To get direct interpolation for Bezier
//https://stackoverflow.com/questions/15505392/implementing-ease-in-update-loop/15506642#15506642
// t is :
// 0.25 = 0.6*t(1-t)^2 + 0.5*t^2(1-t) + t^3
// 0.25 = 0.6 * t *(1 - t) ^ 2 + 0.5 * t^2 * (1 - t) + t ^ 3

// Or
// https://stackoverflow-com.translate.goog/questions/5883264/interpolating-values-between-interval-interpolation-as-per-bezier-curve?_x_tr_sl=en&_x_tr_tl=fr&_x_tr_hl=fr&_x_tr_pto=sc


class cubic_bezier_curve : public bezier_curve_base
{
public:
    cubic_bezier_curve(control_point _cp1, control_point _cp2)
        : _control_point1(_cp1)
        , _control_point2(_cp2)
    {}

    void on_init() override
    {
        c_x = ((-1) * 0 + 3 * _control_point1.x - 3 * _control_point2.x + 1);
        b_x = (3 * 0 - 6 * _control_point1.x + 3 * _control_point2.x);
        a_x = ((-3) * 0 + 3 * _control_point1.x);
        c_y = ((-1) * y_start + 3 * _control_point1.y - 3 * _control_point2.y + y_destination);
        b_y = (3 * y_start - 6 * _control_point1.y + 3 * _control_point2.y);
        a_y = ((-3) * y_start + 3 * _control_point1.y);
    }

    /*
    double process(double x) override
    {
        return y_from_t(calculate_t(x));
    }
    */
private:

    double y_from_t(double t)
    {
        return cubed(1 - t) * y_start
        + 3 * t * squared(1 - t) * _control_point1.y
        + 3 * squared(t) * (1 - t) * _control_point2.y
        + cubed(t) * y_destination;
    }

    double calculate_t(double x)
    {
        if(x == 0) return 0;
        if(x == 1) return 1;
        // Here 0 and 1 stand as P0.x and P3.x
        const double a = -0.0 + 3.0 * _control_point1.x - 3 * _control_point2.x + 1.0;
        const double b = 3.0 * 0.0 - 6 * _control_point1.x + 3.0 * _control_point2.x;
        const double c = -3.0 * 0.0 + 3.0 * _control_point1.x;
        const double d = 0.0 - x;
        return solve_cubic(a, b, c, d);
    }

    double solve_cubic(double a, double b, double c, double d)
    {
        if (a == 0) return solve_quadratic(b, c, d);
        if (d == 0) return 0;

        b /= a;
        c /= a;
        d /= a;
        double q = (3.0 * c - squared(b)) / 9.0;
        double r = (-27.0 * d + b * (9.0 * c - 2.0 * squared(b))) / 54.0;
        double disc = cubed(q) + squared(r);
        double term1 = b / 3.0;

        if (disc > 0) {
            double s = r + std::sqrt(disc);
            s = (s < 0) ? -cubic_root(-s) : cubic_root(s);
            double t = r - std::sqrt(disc);
            t = (t < 0) ? -cubic_root(-t) : cubic_root(t);

            double result = -term1 + s + t;
            if (result >= 0 && result <= 1) return result;
        } else if (disc == 0) {
            double r13 = (r < 0) ? -cubic_root(-r) : cubic_root(r);

            double result = -term1 + 2.0 * r13;
            if (result >= 0 && result <= 1) return result;

            result = -(r13 + term1);
            if (result >= 0 && result <= 1) return result;
        } else {
            q = -q;
            double dum1 = q * q * q;
            dum1 = std::acos(r / std::sqrt(dum1));
            double r13 = 2.0 * std::sqrt(q);

            double result = -term1 + r13 * std::cos(dum1 / 3.0);
            if (result >= 0 && result <= 1) return result;

            result = -term1 + r13 * std::cos((dum1 + 2.0 * M_PI) / 3.0);
            if (result >= 0 && result <= 1) return result;

            result = -term1 + r13 * std::cos((dum1 + 4.0 * M_PI) / 3.0);
            if (result >= 0 && result <= 1) return result;
        }

        throw(std::runtime_error("no result for cubic solve"));
        return 0;
    }

    double solve_quadratic(double a, double b, double c)
    {
        double result = (-b + std::sqrt(squared(b) - 4 * a * c)) / (2 * a);
        if (result >= 0 && result <= 1) return result;

        result = (-b - std::sqrt(squared(b) - 4 * a * c)) / (2 * a);
        if (result >= 0 && result <= 1) return result;

        throw(std::runtime_error("no result for quadratic solve"));
        return 0;
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

    inline virtual double process_all(size_t size, memory_vector<double>::iterator &it) override
    {
        std::vector<double>& res = spl.interpolate_from_points(_control_points, size, point{1.0, 1.0});
        double max = 0.0;
        for(size_t i = 0; i < size; i++)
        {
            if( std::abs(res[i]) > max ) max = std::abs(res[i]);
            if(inverted) res[i] = process_invert(frac(i, size), res[i]);
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
    catmull_rom_spline_curve(point p0, point p3)
        : _cp0(p0)
        , _cp3(p3)
    {}


    void on_init() override
    {
        _cp1 = point(0, y_start);
        _cp2 = point(1, y_destination);
    }


    inline virtual double process_all(size_t size, memory_vector<double>::iterator &it) override
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
            if(inverted) linear_interp = process_invert(x, linear_interp);
            *it = linear_interp;
            ++it;
        }
        return  max;
    }
private:

    std::pair<double, double> process_catmul_rom(double x)
    {
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

    control_point _cp0, _cp3, _cp1, _cp2;
};


}
#endif // CURVE_LIB_H

