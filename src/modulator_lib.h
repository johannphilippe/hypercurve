/*=============================================================================
   Copyright (c) 2022-2024 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#ifndef MODULATOR_LIB_H
#define MODULATOR_LIB_H

#include"curve_lib.h"
#include<cmath>
#include<variant>
#include<type_traits>

namespace hypercurve {
////////////////////////////////////////////////////
// Interpolator is a purely algorithmic curve
// it is used to scale modulators
////////////////////////////////////////////////////
class interpolator
{
public:
    interpolator(double y_start, std::vector<point> pts)
        : itps(std::move(pts))
    {
        // add point 0
        itps.insert(itps.begin(), point(0, y_start));
        for(size_t i = 1; i < itps.size(); ++i)
            itps[i].x += itps[i-1].x;
    }

    virtual double interpolate(double x)
    {
        int seg = which_segment(x);
        double relative_x = relative_position(itps[seg].x, itps[seg + 1].x, x);
        return process(itps[seg].y, itps[seg + 1].y, relative_x);
    }

    virtual double process(double y1, double y2, double x)
    {
        return linear_interpolation(y1, y2, x);
    }

private:
    int which_segment(double x)
    {
        for(size_t i = 0; i < itps.size() - 1; ++i)
        {
            if(x >= itps[i].x && x <  itps[i + 1].x)
                return i;
        }
        return -1;
    }
    std::vector<point> itps;
};

using linear_interpolator = interpolator;

class cubic_interpolator : public interpolator
{
public:
    cubic_interpolator(double y_start ,std::vector<point> pts)
        : interpolator( y_start, std::move(pts) )
    {}
    double process(double y1, double y2, double x) override
    {
        return cubic_interpolation(y1, y2, x);
    }
};

////////////////////////////////////////////////////
// Amplitude can be static (double between 0 and 1)
// or dynamic (interpolator)
////////////////////////////////////////////////////

class amplitude {
public:
    virtual double get_amplitude(double x) {return x;}
};

class amplitude_fixed : amplitude
{
public:
    amplitude_fixed(double d)
        : _amplitude(d)
    {}
    double get_amplitude(double) override {return _amplitude;}
    double _amplitude;
};

class amplitude_interpolated : amplitude
{
public:
    amplitude_interpolated(std::shared_ptr<interpolator> itp)
        : _amplitude(std::move(itp))
    {}
    double get_amplitude(double x) override {return _amplitude->interpolate(x);}
    std::shared_ptr<interpolator> _amplitude;
};


////////////////////////////////////////////
// Modulators are like curve_base
// They are special kind of "curve" that take
// an amplitude parameter (e.g. noise, oscillator)
///////////////////////////////////////////

template<typename Amp>
class modulator_base : public curve_base, public Amp
{
public:
    modulator_base(double amp)
        : Amp(amp)
    {}

    modulator_base(shared_ptr<interpolator> itp)
        : Amp(itp)
    {}

    double process_all(size_t size, memory_vector<double>::iterator &it) override
    {
        double max = 0.0;
        for(size_t i = 0; i < size; ++i)
        {
            const double x = hypercurve::fraction(i, size);
            double res = Amp::get_amplitude(x) * process(x);
            if( std::abs(res) > max) max = std::abs(res);
            *it = res;
            ++it;
        }
        return max;
    }

protected:
};

template<typename Amp = amplitude_fixed>
class noise_modulator : public modulator_base<Amp>
{
public:
    noise_modulator(double amp, size_t precision = 256)
        : modulator_base<Amp>(amp)
        , _precision(precision)
    {}

    noise_modulator(shared_ptr<interpolator> itp, size_t precision = 256)
        : modulator_base<Amp>(itp)
        , _precision(precision)
    {}

    inline double process(double) override
    {
        return ( (rand() % (_precision*2) ) - _precision) / double(_precision);
    }

protected:
    size_t _precision;
};

template<typename Amp>
class sine_modulator : public modulator_base<Amp>
{
public:
    sine_modulator(double amp, double freq)
        : modulator_base<Amp>(amp)
        , _freq(freq)
    {}
    sine_modulator(std::shared_ptr<interpolator> itp, double freq)
        : modulator_base<Amp>(itp)
        , _freq(freq)
    {}

    inline double process(double x) override
    {
        return std::sin(x * M_PI * 2 * _freq);
    }
private:
    double _freq;

};

//////////////////////////////////////////////////
// Chebyshev
// T=1 is stable (basically scaled between -1 and 1)
// T=2 is not, and shouldn't be used if you don't
// know what you're doing
//////////////////////////////////////////////////

template<typename Amp, int T = 1>
class chebyshev_modulator : public modulator_base<Amp>
{
public:
    chebyshev_modulator(double amp, int n_)
        : modulator_base<Amp> (amp)
        , n(n_)
    {}
    chebyshev_modulator(shared_ptr<interpolator> tip, int n_)
        : modulator_base<Amp> (tip)
        , n(n_)
    {}

    inline double process(double x) override
    {

        const double t = std::acos(( x * 2.0) - 1.0 );
        if constexpr(T == 1)
        {
            return std::cos(n * t);
        } else // T == 2
        {
            return std::sin( (n + 1) * t) / std::sin(t);
        }
    }
private:

    const double n;
};


// Be careful, it is cutoff frequencies. The band 1, ends in cutoff 2. So the index 1 corresponds to start index 0, and width 1.
static constexpr int bark_bands_cutoff_frequencies[] = {
    20, 100,200,300,400,510,630, 770, 920, 1080, 1270, 1480, 1720, 2000, 2320, 2700, 3150, 3700, 4400, 5300, 6400, 7700, 9500, 12000, 15500
};

static constexpr int bark_band_widths[] = {
  -1, 80, 100, 100, 100, 110, 120, 140, 150, 160, 190, 210, 240, 280, 320, 380, 450, 550, 700, 900, 1100, 1300, 1800, 2500, 3500
};

static constexpr double bark_width = 15500 - 20;
static constexpr int bark_nbands = sizeof(bark_band_widths) / sizeof(int);
inline void scale_to_bark(double *array, size_t size)
{
    for(size_t i = 0; i < size; i++) {
        double value = array[i] * bark_width + 20;
        for(size_t b = 0; b < bark_nbands; ++b ) {
            if(b+1 > bark_nbands) break;
            if(value >= bark_bands_cutoff_frequencies[b] && value < bark_bands_cutoff_frequencies[b+1]) {
                value *= bark_band_widths[b+1];
                array[i] = value;
                break;
            }
        }
    }
}



}

#endif // MODULATOR_LIB_H
