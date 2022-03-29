#ifndef CURVE_LIB_H
#define CURVE_LIB_H

#include<math.h>

namespace  hypercurve {

class curve_base
{
public:
    curve_base() {}
    virtual ~curve_base() {}

    inline virtual double process(double x) {return x;};

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

}
#endif // CURVE_LIB_H
