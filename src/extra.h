/*=============================================================================
   Copyright (c) 2022 Johann Philippe
   Distributed under the MIT License (https://opensource.org/licenses/MIT)
=============================================================================*/

#pragma once
#ifndef EXTRA_H
#define EXTRA_H

#include<iostream>
#include<unordered_map>
#include"curve_lib.h"
#include"core.h"
#include<memory>

namespace hypercurve {

//////////////////////////////////////////////////////////////////////
// Random and algorithmic curve composition
//////////////////////////////////////////////////////////////////////

// A map of curve types
enum curve_base_index {
    linear_i = 0,
    diocles_i = 1,
    cubic_i = 2,
    power_i = 3,
    hanning_i = 4,
    hamming_i = 5,
    blackman_i = 6,
    gaussian_i = 7,
    toxoid_i = 8,
    catenary_i = 9,
    tightrope_walker_i = 10,
    cubic_bezier_i = 11,
    quadratic_bezier_i = 12,
    cubic_spline_i = 13,
    catmull_rom_spline_i = 14,
    polynomial_i = 15,
    user_defined_i = 16,
    typed_i = 17,
    mouth_i = 18,
    bicorn_i = 19,
    lagrange_polynomial_i = 20,

    // keep that last to get size
    size_i
};

const char *get_curve_base_name (const curve_base_index b)
{
    switch(b) {
    case linear_i: return "linear";
    case diocles_i: return "diocles";
    case cubic_i: return "cubic";
    case power_i: return "power";
    case hanning_i: return "hanning";
    case hamming_i: return "hamming";
    case blackman_i: return "blackman";
    case gaussian_i: return "gaussian";
    case toxoid_i: return "toxoid";
    case catenary_i: return "catenary";
    case tightrope_walker_i: return "tightrope_walker";
    case cubic_bezier_i: return "cubic_bezier";
    case quadratic_bezier_i: return "quadratic_bezier";
    case cubic_spline_i: return "cubic_spline";
    case catmull_rom_spline_i: return "catmull_rom_spline";
    case polynomial_i: return "polynomial";
    case user_defined_i: return "user_defined";
    case typed_i: return "typed";
    case mouth_i: return "mouse";
    case bicorn_i: return "bicorn";
    case lagrange_polynomial_i: return "lagrange_polynomial";
    default: return "";
    }
    return "";
}

std::shared_ptr<curve_base> get_curve_from_index(curve_base_index n,
                                                 std::vector<double> args,
                                                 std::vector<control_point> cps)
{
    switch (n) {
    case linear_i: return share(linear_curve());
    case diocles_i: return share(diocles_curve(args[0]));
    case cubic_i: return share(cubic_curve());
    case power_i: return share(power_curve(args[0]));
    case hanning_i: return share(hanning_curve());
    case hamming_i: return share(hamming_curve());
    case blackman_i: return share(blackman_curve());
    case gaussian_i: return share(gaussian_curve(args[0], args[1]));
    case toxoid_i: return share(toxoid_curve(args[0]));
    case catenary_i: return share(catenary_curve(args[0]));
    case tightrope_walker_i: return share(tightrope_walker_curve(args[0], args[1]));
    case cubic_bezier_i: return share(cubic_bezier_curve(cps[0], cps[1]));
    case quadratic_bezier_i: return share(quadratic_bezier_curve(cps[0]));
    case catmull_rom_spline_i: return share(catmull_rom_spline_curve(cps[0], cps[1]));
    case polynomial_i: return share(polynomial_curve(args));
    case cubic_spline_i: return share(cubic_spline_curve(cps));
    case typed_i: return share(typed_curve(args[0]));
    case mouth_i: return share(mouth_curve());
    case bicorn_i: return share(bicorn_curve(args[0] > 0 ));
    case lagrange_polynomial_i: return share(lagrange_polynomial_curve(cps));
    default: return share(linear_curve());
    }
}

std::pair<std::vector<double>, std::vector<control_point>> random_args_generator(curve_base_index n)
{
    std::vector<double> args;
    std::vector<control_point> cps;

    switch (n) {
    case linear_i: return {{},{}};
    case diocles_i: return {{random<double>(0.50001, 10)}, {}};
    case cubic_i: return {{},{}};
    case power_i: return {{(double)random<int>(1, 32)}, {}};
    case hanning_i: return {{},{}};
    case hamming_i: return {{},{}};
    case blackman_i: return {{},{}};
    case gaussian_i: return {{random<double>(0.1, 4), random<double>(0.00001, 4)}, {}};
    case toxoid_i: return {{random<double>(0.0001, 10)}, {}};
    case catenary_i: return {{random<double>(0.0001, 1.9999)},{}};
    case tightrope_walker_i: {
        double a = random<double>(1, 3);
        double b = a - random<double>(0.01, 0.99);
        return {{a, b}, {}};

    };
    case cubic_bezier_i: {
        std::vector<control_point> cps{
          control_point(random<double>(0, 1), random<double>(0, 1)),
          control_point(random<double>(0, 1), random<double>(0, 1))
        };
        std::sort(cps.begin(),cps.end(), [&](control_point f1,control_point f2) {
            return f1.x < f2.x;
        });

        return {{}, cps};
    };
    case quadratic_bezier_i: return {{}, {control_point(random<double>(0, 1),
                                                        random<double>(0, 1))}};
    case catmull_rom_spline_i: {
        std::vector<control_point> cps {
          control_point(random<double>(0.01, 3) * -1.0, random<double>(0.01, 34) * -1),
          control_point(random<double>(1.01, 3) , random<double>(1.01, 3) ),
        };
        return {{}, cps};
    };
    case polynomial_i: {
        size_t nargs = (rand() % 10) + 1;
        std::vector<double> args(nargs);
        for(size_t i = 0; i < nargs; ++i)
            args[i] = random<double>(0, 10);
        return {args, {}};
        };
    case lagrange_polynomial_i:
    case cubic_spline_i: {
        size_t nargs = (rand() % 4) + 3;
        std::vector<control_point> cps(nargs);
        for(size_t i = 0; i < nargs; ++i)
            cps[i] = control_point(random<double>(0, 1), random<double>(0, 1));
        std::sort(cps.begin(), cps.end(), [&](control_point cp1, control_point cp2) {
           return cp1.x < cp2.x;
        });
        return {{}, cps};
    };
    case typed_i: return {{random<double>(-10, 10)}, {}};
    case mouth_i: return {{},{}};
    case bicorn_i: return {{(double)random<int>(-1, 1)}, {}};

    default: return {{},{}};
    }
}

std::pair<curve, std::string> random_curve_composer( size_t max_segs = 16, int min = 0, int max = 1,
                             size_t definition = 4096, bool envelop = false, bool waveform = false,
                             bool force_curve_type = false, curve_base_index forced = linear_i)
{
    auto gen_curve = [&](){
        return force_curve_type ? forced : static_cast<curve_base_index>(rand() % static_cast<int>(size_i));
    };

    std::string cnames = "";
    size_t nsegs = 1 + (rand() % max_segs);
    if(nsegs < 2) nsegs = 2;


    std::vector<segment> segs(nsegs);
    for(size_t i = 0; i < nsegs; ++i) {
        double frac_size = random<double>(0.1, 1);

        double dest = waveform ? random<double>(-1, 1) : random<double>(0, 1);
        curve_base_index index =  gen_curve();

        while( index == user_defined_i || (envelop && ((index == polynomial_i) || (index == lagrange_polynomial_i) || (index == cubic_spline_i))))
            index = gen_curve();

        cnames += std::string(get_curve_base_name(index)) + "_";
        std::pair<std::vector<double>, std::vector<control_point>> args = random_args_generator(index);
        if(i == (nsegs - 1) && (envelop | waveform) )
            dest = 0;
        segs[i] = segment(frac_size, dest, get_curve_from_index(index, args.first, args.second) );
    }

    double y_start = (envelop | waveform) ? 0 : random<double>(0, 1);
    curve c(definition, y_start, segs);
    c.normalize_y(min, max);
    return {c, cnames};
}

}


#endif // EXTRA_H
