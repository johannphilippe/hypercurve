#ifndef CORE_H
#define CORE_H

#include<math.h>
#include<vector>
#include<memory>
#include"curve_lib.h"
#include"asciiplot/asciiplotter.h"
#include<iostream>

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

    segment(double frac, double y_dest)
        : fractional_size(frac)
        , y_destination(y_dest)
    {}

    virtual double process(std::vector<double>::iterator &it, size_t size)
    {
        return _curve->process_all(size, it);
    }

    virtual void set_y_start(double y)
    {
        y_start = y;
        _curve->init(y_start, y_destination);
    }

    void rescale_x(double factor)
    {
        fractional_size *= factor;
    }

    double fractional_size;
    double y_destination;
    double y_start = 0;
protected:
    std::shared_ptr<curve_base> _curve;
};

/////////////////////////////////////////////////////
// Curve class
////////////////////////////////////////////////////

class curve
{
public:
    curve(size_t definition_, double y_start_, std::vector< std::shared_ptr<segment> > segs_)
        : definition(definition_)
        , y_start(y_start_)
        , segs(segs_)
        , samples(definition)
    {
        init();
        check_total_size();
        process();
    }

    void init()
    {
        segs[0]->set_y_start(y_start);
        for(size_t i = 1; i < segs.size(); i++)
        {
            segs[i]->set_y_start(segs[i - 1]->y_destination);
        }
    }

    void process()
    {
        std::vector<double>::iterator it = samples.begin();
        //
        for(size_t i = 0; i < segs.size(); i++)
        {
            // For each segment, we must give it a real size (size_t), and an iterator position
            size_t seg_size = std::floor(segs[i]->fractional_size * definition);
            double seg_max = segs[i]->process(it, seg_size);
            if(seg_max > max) max = seg_max;
        }
    }

    void normalize_y()
    {

        for(size_t i = 0; i < samples.size(); i++)
        {
            samples[i] = (samples[i] / max);
        }
    }

    void ascii_display(std::string name, std::string label, char marker)
    {
        AsciiPlotter plot(name, 80, 15);
        plot.addPlot(linspace(definition), samples, label, marker);
        plot.legend();
        plot.show();
    }

    double *get_samples() {return samples.data();}
    size_t get_definition() {return definition;}
private:
    void check_total_size()
    {
        double x = 0;
        for(auto & it : segs)
            x += it->fractional_size;
        if( (x > 1.0) ||  (x < 1.0) )
        {
            this->rescale(x);
        }
    }

    void rescale(double x)
    {
        double factor = (1 + (1 - x));
        for(auto & it : segs) {
            it->rescale_x(factor);
        }
    }

    double max = -1.0;
    size_t definition;
    double y_start;
    std::vector< std::shared_ptr<segment> > segs;
    std::vector<double> samples;
};

}

#endif // CORE_H
