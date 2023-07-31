/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#pragma once
#ifndef CURVE_LIB_H
#define CURVE_LIB_H

#include"utilities.h"
#include<cmath>
#include<iostream>
#include<functional>
#include<vector>
#include"cubic_spline.h"

#include"asciiplot/asciiplotter.h"
#include<complex>
typedef std::complex<double> pnt;

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
    inline virtual double process(size_t i, size_t size) {return process(hypercurve::fraction(i, size));}
    inline virtual double process_all(size_t size, memory_vector<double>::iterator &it)
    {
        memory_vector<double>::iterator begin_ptr = it;
        double max = 0.0;

        for(size_t i = 0; i < size; ++i)
        {
            double x = hypercurve::fraction(i, size);

            double res = scale(x);

            if( std::abs(res) > max) max = std::abs(res);

            if(vinverted) res = process_vinvert(x, res);

            *it = res;
            ++it;
        }

        post_processing(begin_ptr);
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

    void post_processing(memory_vector<double>::iterator it)
    {
        if(hinverted) process_hinvert(it);
        if(mirrored) mirror(it, definition, y_start, y_destination);
    }

    // Override this one insted of init to avoid y_start and y_destination
    // affectation repetition
    inline virtual void on_init() {}

    // Allows inversion of curve (y symetry on a linear x_start/x_end axis).
    bool hinverted = false;
    bool vinverted = false;
    bool mirrored = false;
protected:

    inline virtual double scale(double x)
    {
        double y;
        if(y_start > y_destination)
            y = process(1.0 - x) ;
        else
            y = process(x); // * abs_diff + offset;

        y  = y * abs_diff + offset;
        return y;
    }

    // Vertical inversion
    inline double process_vinvert(double x, double y)
    {
        const double lin = linear_interpolation(y_start, y_destination, x);
        //std::cout << std::fixed << " lin  " << lin << " y " << y << " res " << (lin + (lin - y)) << std::endl;
        return lin + (lin - y) ;
    }

    inline void process_hinvert(memory_vector<double>::iterator it)
    {
        auto endit = it + definition;
        auto itrev = it;
        abs_diff = std::abs(y_destination - y_start);
        std::reverse(it.ptr, endit.ptr);
        for(size_t i = 0; i < definition; ++i)
        {
           *(it.ptr) = (abs_diff + y_start) - (*(it.ptr));
            ++it;
        }
    }

    size_t definition;
    double y_start, y_destination, abs_diff, offset;
};

using linear_curve = curve_base;

// Allows you to make a symetry on a x_start/y_start x_end/y_end linear axis
inline std::shared_ptr<curve_base> vinvert(std::shared_ptr<curve_base> cb)
{
    cb->vinverted = !cb->vinverted;
    return cb;
}
inline std::shared_ptr<curve_base> vsymmetry(std::shared_ptr<curve_base> cb) {return vinvert(cb);}
inline std::shared_ptr<curve_base> vsym(std::shared_ptr<curve_base> cb) {return vinvert(cb);}

inline std::shared_ptr<curve_base> hinvert(std::shared_ptr<curve_base> cb)
{
    cb->hinverted = !cb->hinverted;
    return cb;
}
inline std::shared_ptr<curve_base> hsymmetry(std::shared_ptr<curve_base> cb) {return hinvert(cb);}
inline std::shared_ptr<curve_base> hsym(std::shared_ptr<curve_base> cb) {return hinvert(cb);}

inline std::shared_ptr<curve_base> mirror(std::shared_ptr<curve_base> cb)
{
    cb->mirrored = !cb->mirrored;
    return cb;
}

//////////////////////////////////////////////////
// Curve Library
//////////////////////////////////////////////////



//////////////////////////////////////////////////
// Logarithmic and exponential approximations
//////////////////////////////////////////////////
class logarithmic_curve : public curve_base
{
public:
    logarithmic_curve()
    {
    }
    inline double process(double x) override
    {
        return scaled_log_point<double>(x + offset, offset, 1 + offset) / (1 + offset);
    }

    constexpr static const double offset = 0.1;
};

class exponential_curve : public logarithmic_curve
{
public:
    exponential_curve()
    {
        mirrored = !mirrored;
    }
};



class diocles_curve : public curve_base
{
public:
    diocles_curve(double a_)
        : a(a_)
        , compensation(1. / process_diocles(1.0) )
    {
        if(a_ <= 0.5) throw(std::runtime_error("Diocles curve : 'a' parameter must be > to 0.5"));
    }
    inline double process(double x) override
    {
        return process_diocles(x) * (compensation);
    }

private:

    inline double process_diocles(double x)
    {
         return std::sqrt( std::pow(x, 3) / (2 * a - x) );
    }

    double a, compensation;
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
        return std::pow(x, exponent) ;
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
        return  (hamming(x * definition, definition * 2) - hamming_scaling_constant) * hamming_scaling_factor;
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
class mouth_curve : public curve_base
{
public:
    inline double process(double x) override
    {
        const double _x = x -1;
        return std::sqrt( std::pow( (a*a) - (_x*_x), 3 ) / std::pow(a, 4) );
    }

    constexpr static const double a = 1.0;
};
using kiss_curve = mouth_curve;

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
        std::cout <<"conchal : " << ((c*c) + (a*a) - (x*x) )
                    * ((c*c) - (a*a) + (x*x)) << "   " << (x+a) << std::endl;
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
// The more a is superior to b, the more is is "straight"
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
        , height(toxoid_process( 1))
    {
    }

    double process(double x) override
    {
        return toxoid_process(x) / height;
    }

protected:

    inline double toxoid_process(double x)
    {
        return std::sqrt(x * std::pow(x - a/2, 2));
        //return std::sqrt( ( (x*x) * (x - a/2) ) / a  );
    }

    double a, height;
};
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
        memory_vector<double>::iterator begin_ptr = it;
        double max = 0.0;
        size_t cnt = 0;
        std::pair<double, double> r1, r2;
        r1 = process_bezier(double(cnt) / double(size));
        r2 = process_bezier(double(cnt + 1) / double(size));

        for(size_t i = 0; i < size ; i++)
        {
            double x = hypercurve::fraction(i, size);
            //const double x = double(i) / double(size);
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
            if(vinverted) linear_interp = process_vinvert(x, linear_interp);
            *it = linear_interp;
            ++it;
        }
        post_processing(begin_ptr);
        return max;
    }

    inline virtual double process(double) override
    {
        throw(std::runtime_error("Unimplemented for Bezier curve"));
    }

    // This one should be implemented instead of the above one (if wanted to process single bezier point)
    inline virtual double process(size_t i, size_t size) override
    {
        return process(hypercurve::fraction(i, size));
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

    // Internal allocator for Spline aux mem
    cubic_spline_curve(std::vector<point> cp)
        : _control_points( cp.size() + 2 )
    {
        if(_control_points.size() < 3)
            throw(std::runtime_error("Control point list size must be >= 3"));
        std::move(cp.begin(), cp.end(), _control_points.begin() + 1);
        int n = _control_points.size() - 1;
        _spline_mem_size = (5 * n) + (3 * ( n + 1 )) + (n - 1)  + ( (n - 1) * n);
        _spline_memory.resize(_spline_mem_size);
        _spline_memory_ptr = _spline_memory.data();
    }

    // Internal allocator for Spline aux mem
    cubic_spline_curve(memory_vector<point> cp)
        : _control_points(  cp.size() + 2 )
    {
        if(_control_points.size() < 3)
            throw(std::runtime_error("Control point list size must be >= 3"));

        int n = _control_points.size() - 1;
        _spline_mem_size = (5 * n) + (3 * ( n + 1 )) + (n - 1)  + ( (n - 1) * n);
        _spline_memory.resize(_spline_mem_size);
        _spline_memory_ptr = _spline_memory.data();
        std::move(cp.begin(), cp.end(), _control_points.begin() + 1);
    }


    // External allocator constructor
    cubic_spline_curve(point *pts, size_t pts_size, double *spline_memory, size_t spline_mem_size)
        : _control_points(pts, pts_size)
        , _spline_memory_ptr(spline_memory)
        , _spline_mem_size(spline_mem_size)
    {
        if(pts_size < 3)
            throw(std::runtime_error("Control point list size must be >= 3"));
    }

    ~cubic_spline_curve(){
    }

    void on_init() override
    {

        _control_points[0] = control_point(0, y_start);
        _control_points[_control_points.size() - 1] = control_point(1, y_destination);
    }

    inline virtual double process_all(size_t size, memory_vector<double>::iterator &it) override
    {
        memory_vector<double>::iterator begin_ptr = it;
        memory_vector<double>::iterator end_ptr = it;

        double *data_ptr = it.get();
        spl.process(data_ptr, definition, _control_points, _spline_memory_ptr, _spline_mem_size);
        double max = 0.0;

        for(size_t i = 0; i < size; i++)
        {
            if( std::abs(*it) > max ) max = std::abs(*it);
            if(vinverted) *it = process_vinvert(hypercurve::fraction(i, size), *it);
            ++it;
        }

        post_processing(begin_ptr);
        return max;
    }

private:
    cubic_spline<double> spl;
    memory_vector<point> _control_points;

    // Vector, automatically destroyed when using internal allocator
    std::vector<double> _spline_memory;
    // Ptr to allocated memory, whether it is allocated externally or internally
    double* _spline_memory_ptr = nullptr;
    size_t _spline_mem_size;
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
        memory_vector<double>::iterator begin_ptr = it;
        size_t cnt = 0;
        std::pair<double, double> r1, r2;
        r1 = process_catmul_rom(0);
        r2 = process_catmul_rom(1.0 / double(size));
        double max = 0.0;
        for(size_t i =0; i < size; i++)
        {
            double x = hypercurve::fraction(i, size);

            //const double x = double(i) / double(size);
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
            if(vinverted) linear_interp = process_vinvert(x, linear_interp);
            *it = linear_interp;
            ++it;
        }

        post_processing(begin_ptr);
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


class lagrange_polynomial_curve : public curve_base
{

public:
    lagrange_polynomial_curve(std::vector<control_point> pts)
        : c_pts(pts.size() + 2)
    {
        std::copy(pts.begin(), pts.end(), c_pts.data() + 1);
    }
    lagrange_polynomial_curve(memory_vector<control_point> pts)
        : c_pts(pts.size() + 2)
    {
        std::copy(pts.begin(), pts.end(), c_pts.data() + 1);
    }

    // Used for preallocated situations (e.g. Csound)
    lagrange_polynomial_curve(control_point *pts, size_t size)
        : c_pts(pts, size)
    {}

    void on_init() override
    {
        c_pts[0] = control_point(0, y_start);
        c_pts[c_pts.size() - 1] = control_point(1, y_destination);
    }

    inline double process(double x) override
    {
        return process_lagrange(x);
    }

protected:
    double process_lagrange(double x)
    {
        yp = 0;
        for(size_t i = 0; i < c_pts.size(); ++i)
        {
            p = 1;
            for(size_t j = 0; j < c_pts.size(); ++j)
            {
                if(j != i)
                    p = p * (x - c_pts[j].x)/( c_pts[i].x  - c_pts[j].x);
            }
            yp = yp + p * c_pts[i].y;
        }
        return yp;
    }

    double p, yp;
    memory_vector<control_point> c_pts;
};

/**
  ***************************************************************
 * Easing curves from easings.net
 * - Missing Quad, Cubic, Quart, Quint, Expo
 ****************************************************************
*/


class ease_in_sine : public curve_base
{
public:
    inline double process(double x) override
    {
        return 1.0 - std::cos((x * M_PI) / 2);
    }
};

class ease_out_sine : public curve_base
{
public:
    inline double process(double x) override
    {
        return std::sin((x * M_PI) / 2);
    }
};

class ease_inout_sine : public curve_base
{
public:
    inline double process(double x) override
    {
        return -(std::cos(M_PI * x) - 1.0) / 2.0;
    }
};

class ease_in_circ : public curve_base
{
public:
    inline double process(double x) override
    {
        return 1 - std::sqrt(1 - std::pow(x, 2));
    }
};

class ease_out_circ : public curve_base
{
public:
    inline double process(double x) override
    {
        return std::sqrt(1 - std::pow(x - 1, 2));
    }
};

class ease_inout_circ : public curve_base
{
public:
    inline double process(double x) override
    {
        return x < 0.5
          ? (1 - std::sqrt(1 - std::pow(2 * x, 2))) / 2
          : (std::sqrt(1 - std::pow(-2 * x + 2, 2)) + 1) / 2;
    }
};


class ease_in_back : public curve_base
{
public:
    inline double process(double x) override
    {
        return c3 * x * x * x - c1 * x * x;
    }
    constexpr static const double c1 = 1.70158;
    constexpr static const double c3 = c1 + 1;
};
class ease_out_back : public curve_base
{
public:
    inline double process(double x) override
    {
        return 1 + c3 * std::pow(x - 1, 3) + c1 * std::pow(x - 1, 2);
    }

    constexpr static const double c1 = 1.70158;
    constexpr static const double c3 = c1 + 1;
};
class ease_inout_back : public curve_base
{
public:
    inline double process(double x) override
    {
        return x < 0.5
          ? (std::pow(2 * x, 2) * ((c2 + 1) * 2 * x - c2)) / 2
          : (std::pow(2 * x - 2, 2) * ((c2 + 1) * (x * 2 - 2) + c2) + 2) / 2;
    }

    constexpr static const double c1 = 1.70158;
    constexpr static const double c2 = c1 * 1.525;

};

class ease_in_elastic : public curve_base
{
public:
    inline double process(double x) override
    {
        return x == 0
          ? 0
          : x == 1
          ? 1
          : -std::pow(2, 10 * x - 10) * std::sin((x * 10 - 10.75) * c4);
    }

    constexpr static const double c4 = (2 * M_PI) / 3;
};

class ease_out_elastic : public curve_base
{
public:
    inline double process(double x) override
    {

        return x == 0
          ? 0
          : x == 1
          ? 1
          : std::pow(2, -10 * x) * std::sin((x * 10 - 0.75) * c4) + 1;
    }

    constexpr static const double c4 = (2 * M_PI) / 3;
};

class ease_inout_elastic : public curve_base
{
public:
    inline double process(double x) override
    {
        return x == 0.0
          ? 0.0
          : x == 1.0
          ? 1.0
          : x < 0.5
          ? -(std::pow(2, 20 * x - 10) * std::sin((20 * x - 11.125) * c5)) / 2.0
          : (std::pow(2, -20 * x + 10) * std::sin((20 * x - 11.125) * c5)) / 2.0 + 1.0;
    }

    constexpr static const double c5 = (2 * M_PI) / 4.5;
};



class ease_out_bounce : public curve_base
{
public:

    inline double process(double x) override
    {
        return process_value(x);
    }

    static double process_value(double x ) {
        if (x < 1.0 / d1) {
            return n1 * x * x;
        } else if (x < 2 / d1) {
            x -= (1.5 / d1);
            return n1 * (x) * x + 0.75;
        } else if (x < 2.5 / d1) {
            x -= 2.25 / d1;
            return n1 * (x) * x + 0.9375;
        } else {
            x -= (2.625 / d1);
            return n1 * (x) * x + 0.984375;
        }
    }
    constexpr static const double n1 = 7.5625;
    constexpr static const double d1 = 2.75;
};


class ease_in_bounce : public curve_base
{
public:
    inline double process(double x) override
    {
        return 1.0 - ease_out_bounce::process_value(1.0 - x);
    }
};

class ease_inout_bounce : public curve_base
{
public:
    inline double process(double x) override
    {
        return x < 0.5
          ? (1.0 - ease_out_bounce::process_value(1.0 - 2.0 * x)) / 2.0
          : (1.0 + ease_out_bounce::process_value(2.0 * x - 1)) / 2.0;
    }
};


}
#endif // CURVE_LIB_H

