/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/
#pragma once
#ifndef CORE_H
#define CORE_H

#include<cmath>
#include<vector>
#include<memory>
#include"curve_lib.h"
#include"asciiplot/asciiplotter.h"
#include<iostream>
#include<utility>
namespace hypercurve {
///////////////////////////////////////////////////
// The segment class
////////////////////////////////////////////////////

class segment
{
public:
    segment(double frac, double y_dest, std::shared_ptr<curve_base> c)
        : fractional_size(frac)
        , y_destination(y_dest)
        , _curve( std::move(c) )
    {}

    segment() {}

    virtual double process(memory_vector<double>::iterator &it, size_t size)
    {
        return _curve->process_all(size, it);
    }

    virtual void init(double y_start_, size_t definition)
    {
        y_start = y_start_;
        _curve->init(y_start, y_destination, definition);
    }
    void rescale_x(double factor)
    {
        fractional_size *= factor;
    }

    double fractional_size;
    double y_destination;
    double y_start = 0;
    std::shared_ptr<curve_base> _curve;
protected:
};

/////////////////////////////////////////////////////
// Curve class
////////////////////////////////////////////////////

class curve
{
public:
    curve(size_t definition_, double y_start_, std::vector< segment > segs_)
        : definition(definition_)
        , y_start(y_start_)
        , segs( std::move(segs_) )
        , samples(definition)
    {
        check_total_size();
        init();
        process();
    }

    curve() {}
    virtual  ~curve() {}

    virtual void init()
    {
        segs[0].init(y_start, std::round(segs[0].fractional_size * definition) );
        for(size_t i = 1; i < segs.size(); i++)
        {
            segs[i].init(segs[i - 1].y_destination, std::round(segs[i].fractional_size * definition) );
        }
    }

    void process()
    {
       memory_vector<double>::iterator it = samples.begin();
        for(size_t i = 0; i < segs.size(); i++)
        {
            // For each segment, we must give it a real size (size_t), and an iterator position
            size_t seg_size = std::round(segs[i].fractional_size * definition);
            double seg_max = segs[i].process(it, seg_size);
            if(std::abs(seg_max)  > max) max = std::abs(seg_max);
        }
        find_extremeness();
    }

    void normalize_y(double target_min, double target_max)
    {
        find_extremeness();

        for(size_t i = 0; i < samples.size(); i++)
        {
            samples[i] = ((samples[i] - min ) / ambitus )  * std::abs(target_max - target_min) + target_min;
        }
    }

    std::pair<double, double> find_extremeness()
    {
        min = samples[0], max = samples[0];
        for(auto & it : samples)
        {
            if(it < min)
                min = it;
            if(it > max)
                max = it;
        }
        ambitus = std::abs(max - min);
        return std::make_pair(min, max);
    }

    void ascii_display(std::string name, std::string label, char marker)
    {
        AsciiPlotter plot(name, 80, 15);
        plot.addPlot(std::vector<double>(samples.data(), samples.data() + definition), label, marker);
        plot.legend();
        plot.show();
    }

    double *get_samples() {return samples.data();}
    double get_sample_at(size_t i) {return samples[i];}
    size_t get_definition() {return definition;}

    // Operators

    double& operator[](size_t index)  { return samples[index];}

    curve& operator *=(curve &other)
    {
        for(size_t i = 0; i < samples.size(); ++i)
            samples[i] *= other.get_sample_at(i);
        return *this;
    }
    curve &operator *=(double k)
    {
        for(auto & it : samples)
            it *= k;
        return *this;
    }

    curve& operator /=(curve &other)
    {
        for(size_t i = 0; i < samples.size(); ++i)
            samples[i] /= other.get_sample_at(i);
        return *this;
    }
    curve &operator /=(double k)
    {
        for(auto & it : samples)
            it /= k;
        return *this;
    }

    curve operator *(curve &c2)
    {
        curve c(*this);
        c *= c2;
        return c;
    }

    curve operator *(double k)
    {
        curve c(*this);
        c *= k;
        return c;
    }

    curve operator /(curve &c2)
    {
        curve c(*this);
        c /= c2;
        return c;
    }

    curve operator / (double k)
    {
        curve c(*this);
        c /= k;
        return c;
    }

    curve& operator+=(double k)
    {
        for(auto & it : samples)
            it += k;
        return *this;
    }

    curve& operator +=(curve &other)
    {
        for(size_t i = 0; i < samples.size(); ++i)
            samples[i] += other.get_sample_at(i);
        return *this;
    }

    curve operator +(double k)
    {
        curve c(*this);
        c += k;
        return c;
    }
    curve operator +(curve &other)
    {
        curve c(*this);
        c += other;
        return c;
    }

    curve& operator -=(double k)
    {
        for(auto & it : samples)
            it -= k;
        return *this;
    }
    curve& operator -=(curve &other)
    {
        for(size_t i = 0; i < samples.size(); ++i)
            samples[i] -= other.get_sample_at(i);
        return *this;
    }

    curve operator -(double k)
    {
        curve c(*this);
        c -= k;
        return c;
    }

    curve operator -(curve &other)
    {
        curve c(*this);
        c -= other;
        return c;
    }

protected:
    void check_total_size()
    {
        double x = 0;
        for(auto & it : segs)
            x += it.fractional_size;
        if( (x > 1.0) ||  (x < 1.0) )
        {
            this->rescale(x);
        }
    }

    void rescale(double x)
    {
        double factor = (1. / x);
        for(auto & it : segs) {
            it.rescale_x(factor);
        }
    }

    double max = 0.0, min = 0.0, ambitus = 0.0;
    size_t definition;
    double y_start;
    memory_vector< segment > segs;
    memory_vector<double> samples;
};

}
#endif // CORE_H
