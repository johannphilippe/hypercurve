#include <iostream>
#include<sndfile.hh>
#include<math.h>
#include<vector>
#include<memory>
using namespace std;


/*
 *  Call curve with full size
 *
 * 		// definition 	// starty
 * curve(16384, 		0.1,
 * {		// size		//toy	// mode (default to linear)
 * 	segment(4/8,	, 0.5	 	cubic()),
 * });
*/

class curve_base
{
public:
    curve_base() {}
    virtual ~curve_base() {}

    virtual double process(double x) = 0;

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

double pos(double x)
{
    return (x > 0) ? x : 0;
}

class segment
{
public:
    segment(double frac, double y_dest, std::shared_ptr<curve_base> c)
        : fractional_size(frac)
        , y_destination(y_dest)
        , _curve(c )
    {}

    void process(std::vector<double>::iterator it, size_t size, double y_start)
    {
        // ambitus
        const double y_diff = std::abs(y_destination - y_start);
        const double inv_diff = 1 - y_diff;
        std::cout << "diff : " << y_diff << " " << inv_diff << std::endl;

        // diff btw start of this curve and its max point
        // diff btw ystart and scale
        for(size_t i = 0; i < size; i++)
        {
            // Get ymin, ymax to scale.
            double res = _curve->process(double(i) / double(size)) ;
            if(y_start > y_destination) res = 1 - res;
            res *= y_diff;
            res += inv_diff;
            *it = res;
            it++;
        }
    }

    const double fractional_size;
    const double y_destination;
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
            std::cout << "seg iterate " << i << std::endl;
            // For each segment, we must give it a real size (size_t), and an iterator position
            size_t seg_size = std::floor(segs[i].fractional_size * definition);
            std::cout << "y_start : " << y_from << std::endl;
            segs[i].process(it, seg_size, y_from);
            y_from = segs[i].y_destination;
            it+= seg_size;
        }
    }

    void normalize_y()
    {

    }

    double *get_samples() {return samples.data();}
private:
    // TODO : implement an X rescale mode to equal zero.

        // Check if sum of all seg size is equal to 0
    void check_total_size()
    {
        double x = 0;
        for(auto & it : segs)
            x += it.fractional_size;
        if(x == 1.0)
        {
            // ALL ok

        } else if(x > 1.0)
        {
            // Should raise error or rescale everything.
            throw(std::runtime_error("Curve components goes across 1.0 fractional. Use rescale mode to rescale"));
        } else if(x < 1.0)
        {
            // Add a linear segment to 0 or rescale
        }
    }


    size_t definition;
    double y_start;
    std::vector<segment> segs;
    std::vector<double> samples;

};

double frac(double x1, double x2)
{
    return double(x1) / double(x2);
}

template<typename T>
std::shared_ptr<T> share(T t)
{
    return std::make_shared<T>(t);
}

int main()
{
    SndfileHandle sndin("/home/johann/Documents/tmp/sweep.wav");

    SndfileHandle sndout("/home/johann/Documents/csound/diocles.wav", SFM_WRITE, sndin.format(), 1, 48000);

    int def = 16384;
    curve c(def, 0
            , 	{
                segment(frac(1,4), 1.0, share(diocles_curve(1))),
                segment(frac(1,4), 0.5, share(cubic_curve())),
                segment(frac(1,4), 1, share(diocles_curve(1))),
                segment(frac(1,4), 0, share(diocles_curve(1)))
                }
            );

    curve c2(def, 0, {segment(frac(1,1), 1.0, share(diocles_curve(1)))});

    sndout.writef(c.get_samples(), def);

    return 0;
}
