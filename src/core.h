#ifndef CORE_H
#define CORE_H

#include<math.h>
#include<vector>
#include<memory>
#include"curve_lib.h"
#include<iostream>

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

class segment
{
public:
    segment(double frac, double y_dest, std::shared_ptr<curve_base> c)
        : fractional_size(frac)
        , y_destination(y_dest)
        , _curve( c )
    {}

    double process(std::vector<double>::iterator it, size_t size, double y_start)
    {
        // ambitus
        const double y_diff = std::abs(y_destination - y_start);
        const double inv_diff = std::min(y_destination, y_start);// 1 - y_diff;

        double max = -1;
        // diff btw start of this curve and its max point
        // diff btw ystart and scale
        for(size_t i = 0; i < size; i++)
        {
            // Get ymin, ymax to scale.
            double res = _curve->process(double(i) / double(size)) ;
            // turn it to get the right direction
            if(y_start > y_destination) res = 1 - res;
            // scale it 		then add offset
            res = (res * y_diff) + inv_diff;
            if(res > max) max = res;
            *it = res;
            it++;
        }
        return max;
    }

    void rescale_x(double factor)
    {
        fractional_size *= factor;
    }

    double fractional_size;
    double y_destination;
private:
    std::shared_ptr<curve_base> _curve;
};


class curve
{
public:
    curve(size_t definition_, double y_start_, std::vector<segment> segs_)
        : definition(definition_)
        , y_start(y_start_)
        , segs(segs_)
        , samples(definition)
    {
        check_total_size();
        process();
    }

    void process()
    {
        double y_from = y_start;

        std::vector<double>::iterator it = samples.begin();
        //
        for(size_t i = 0; i < segs.size(); i++)
        {
            // For each segment, we must give it a real size (size_t), and an iterator position
            size_t seg_size = std::floor(segs[i].fractional_size * definition);
            double seg_max = segs[i].process(it, seg_size, y_from);
            if(seg_max > max) max = seg_max;
            y_from = segs[i].y_destination;
            it+= seg_size;
        }
    }

    void normalize_y()
    {

        for(size_t i = 0; i < samples.size(); i++)
        {
            samples[i] = (samples[i] / max);
        }
    }

    double *get_samples() {return samples.data();}
    size_t get_definition() {return definition;}
private:
    // TODO : implement an X rescale mode to equal zero.

        // Check if sum of all seg size is equal to 0
    void check_total_size()
    {
        double x = 0;
        for(auto & it : segs)
            x += it.fractional_size;
        if( (x > 1.0) ||  (x < 1.0) )
        {
            // Should raise error or rescale everything.
            this->rescale(x);
            //throw(std::runtime_error("Curve components goes across 1.0 fractional. Use rescale mode to rescale"));
        }
    }

    void rescale(double x)
    {
        double factor = (1 + (1 - x));
        for(auto & it : segs) {
            it.rescale_x(factor);
        }
    }

    double max = -1.0;
    size_t definition;
    double y_start;
    std::vector<segment> segs;
    std::vector<double> samples;

};

}

#endif // CORE_H
