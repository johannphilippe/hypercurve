#ifndef CORE_H
#define CORE_H

#include<math.h>
#include<vector>
#include<memory>
#include"curve_lib.h"
#include"asciiplot/asciiplotter.h"
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

    virtual double process(std::vector<double>::iterator it, size_t size)
    {
        // ambitus
        const double y_diff = std::abs(y_destination - y_start);
        const double offset = std::min(y_destination, y_start);// 1 - y_diff;

        double max = -1;
        for(size_t i = 0; i < size; i++)
        {
            // Get ymin, ymax to scale.
            double res = _curve->process(double(i) / double(size)) ;
            // turn it to get the right direction
            if(y_start > y_destination) res = 1 - res;
            // scale it 		then add offset
            res = (res * y_diff) + offset;
            if(res > max) max = res;
            *it = res;
            ++it;
        }
        return max;
    }

    virtual void set_y_start(double y)
    {
        y_start = y;
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

///////////////////////////////////////////////////
// Bezier segment
////////////////////////////////////////////////////
class bezier_segment : public segment
{
public:
    bezier_segment(double frac, double y_dest, std::shared_ptr<bezier_curve_base> crv)
        : segment(frac, y_dest)
        , _curve( std::move(crv) )
    {}

    void set_y_start( double y) override
    {
        y_start = y;
        _curve->init(y_start, y_destination);
    }

    double process(std::vector<double>::iterator it, size_t size) override
    {
        double max = -1;

        size_t cnt = 0;
        std::pair<double, double> r1, r2;
        r1 = _curve->process_bezier(double(cnt) / double(size));
        r2 = _curve->process_bezier(double(cnt + 1) / double(size));

        for(size_t i = 1; i <= size; i++)
        {
            const double x = double(i) / double(size);
            while(cnt < size)
            {
                if( (r1.first <= x ) && (r2.first >= x))
                    break;
                r1 = _curve->process_bezier( double(cnt) / double(size) );
                r2 = _curve->process_bezier( double(cnt + 1) / double(size) );
                ++cnt;
            }

            double relative_x = relative_position(r1.first, r2.first, x);
            double linear_interp = linear_interpolation(r1.second, r2.second, relative_x);
            if(linear_interp > max) max = linear_interp;
            *it = linear_interp;
            ++it;
        }

        return max;
    }

protected:
    std::shared_ptr<bezier_curve_base> _curve;
};

///////////////////////////////////////////////////
// Spline segment
////////////////////////////////////////////////////
class spline_segment : public segment
{
public:
    spline_segment(double frac, double y_dest, std::shared_ptr<spline_curve_base> cp)
        : segment(frac, y_dest)
        , _curve(std::move(cp))
    {
    }

    void set_y_start( double y) override
    {
        y_start = y;
        _curve->init(y_start, y_destination);
    }

    double process(std::vector<double>::iterator it, size_t size) override
    {
        std::vector<double> &res = _curve->interpolate(size);
        double max = -1;
        for(size_t i = 0; i < size; i++)
        {
            if(res[i] > max) max = res[i];
            *it = limit(0, 1, res[i]);
            ++it;
        }
        return max;
    }
private:
    std::shared_ptr<spline_curve_base> _curve;
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
            //throw(std::runtime_error("Curve components goes across 1.0 fractional. Use rescale mode to rescale"));
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
