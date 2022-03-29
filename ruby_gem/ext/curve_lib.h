#ifndef CURVE_LIB_H
#define CURVE_LIB_H

#include<math.h>

namespace  hypercurve {

class curve_base
{
public:
    curve_base() {}
    virtual ~curve_base() {}

    virtual double process(double x) = 0;

};

class linear_curve : public virtual curve_base
{
public:
    double process(double x) override
    {
        return x;
    }
};

class diocles_curve : public virtual curve_base
{
public:
    diocles_curve(const double a_)
        : a(a_)
    {}
    double process(double x) override
    {
        return std::sqrt( std::pow(x, 3) / (2 * a - x) );
    }

private:
    const double a;
};

class cubic_curve: public virtual curve_base
{
public:
    double process(double x) override
    {
        return std::pow(x, 3);
    }
};

}
#endif // CURVE_LIB_H
