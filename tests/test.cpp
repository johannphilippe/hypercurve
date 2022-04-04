#include <iostream>
#include<sndfile.hh>
#include<memory.h>
#include"../src/hypercurve.h"
#include"sndfile.hh"
#include<chrono>
using namespace std;

struct timer
{
    timer()
    {
        t1 = std::chrono::high_resolution_clock::now();
        first = t1;
    }

    void time_since_last(std::string label = "")
    {
        t2 = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count();
        std::cout << label << " : " << double(dur) / 1000000.0  << " ms " << std::endl;
        t1 = t2;
    }

    void time_since_start(std::string label = "")
    {
        t2 = std::chrono::high_resolution_clock::now();
        auto dur = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - first).count();
        std::cout << label << " : " << double(dur) / 1000000.0  << " ms " << std::endl;
    }

    void reset()
    {
        t1 = std::chrono::high_resolution_clock::now();
    }


    std::chrono::high_resolution_clock::time_point first, t1, t2;
};

using namespace hypercurve;

int main()
{
    int def = 16384;

    const int fmt = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
    SndfileHandle sf("test_hypercurveMod.wav", SFM_WRITE, fmt, 1, 48000);
    SndfileHandle sf2("test_hypercurveModulated.wav", SFM_WRITE, fmt, 1, 48000);

    timer t;

    // Simple composite curve
    curve c(def, 0
            , 	{
                segment(frac(1,4), 1.0, share(diocles_curve(1))),
                segment(frac(1,4), 0.5, share(cubic_curve())),
                segment(frac(1,4), 1, share(diocles_curve(1))),
                segment(frac(1,4), 0, share(diocles_curve(1)))
                }
            );

    t.time_since_last("composite curve");

    c.ascii_display("Composite diocles and cubic curve", "y = composite(x)", '*');

    t.reset();

    // One segment curve
    curve c2(def, 0, {segment(frac(1,1), 1.0, share(diocles_curve(1)))});

    t.time_since_last("one segment curve");

    c2.ascii_display("One segment diocles curve", "y = oneseg(x)", '*');
    double div = 4;

    t.reset();
    // Bezier quadratic (one control point)
    curve c3(def, 0, {
                 segment(frac(1,div), 1.0, share(quadratic_bezier_curve({0.1, 0.9}))),
                 segment(frac(1,div), 0.5, share(quadratic_bezier_curve({0.66, 0.1}))),
                 segment(frac(1,div), 0.8, share(quadratic_bezier_curve({0.9, 0.9}))),
                 segment(frac(1, div), 0.1, share(quadratic_bezier_curve({0.5, 0})))
             });

    t.time_since_last("Bezier quadratic 4 segments");

    c3.ascii_display("Bezier quadratic", "y = quadratic_bezier(x)", '*');

    t.reset();
    // Bezier cubic (taking two control points)
    curve c4(def, 0, {
                 segment(frac(1,1), 1	, share(cubic_bezier_curve({0.2, 0.8}, {0.8, 0.2})))
             });

    t.time_since_last("Bezier cubic");

    c4.ascii_display("Bezier Cubic", "y = cubic_bezier(x)", '*');

    t.reset();
    // homemade exponent
    curve c5(def, 0, {
                segment(1, 1, share(exponent_curve(9)))
             });

    t.time_since_last("Homemade exponent");
    c5.ascii_display("Smooth exponent", "y = hsmooth(x)", '*');

    // Cubic spline
    t.reset();
    curve c6(def, 0, {
                 segment(1, 1, {
                     share(cubic_spline_curve({
                        point(0.2, 0.8),
                        point(0.2, 0.2),
                        point(0.3, 0.9),
                        point(0.3, 1)
                     }))
                 })
             });
    t.time_since_last("Cubic spline");
    c6.ascii_display("Cubic spline", "y = cspline(X)", '*');

    // CatmullRom
    t.reset();
    curve c7(def, 0, {
                 segment(frac(1,2), 1, {
                     share(catmull_rom_spline_curve(
                        point(-2, -0.5),
                        point(2, 0.2)))
                 }),
                 segment(frac(1,2), 0, {
                     share(catmull_rom_spline_curve(
                     point(-1, 3),
                     point(3, -5)
                     ))
                 })
             });
    t.time_since_last("Catmull Rom spline");

    c7.ascii_display("CatmullRom", "y = catmulltom(X)", '*');

    t.reset();
    curve	 w(def, 0, {
                  segment(frac(1,3), 1, share(cubic_curve())),
                  segment(frac(1,3), -1, share(cubic_curve())),
                  segment(frac(1,3), 0, share(cubic_curve())),

               });
    t.time_since_last("Waveform");

    w.ascii_display("waveform", "w=cubic(x)", '.');


    curve to_modulate(16384, 0, {
                         segment(frac(1,2), 1, share(diocles_curve(1))),
                         segment(frac(1,2), 0, share(diocles_curve(1))),

                      });

    curve mod(16384, 0, {
               segment(1, 1, share(chebyshev_modulator<amplitude_fixed>(0.1, 15)))
            });

    curve modulated = to_modulate + mod;
    mod.ascii_display("modulator" ,"mod(x)", '*');

    modulated.ascii_display("Modulated", "mod * catmullrom", '*');


    sf.writef(mod.get_samples(), def);
    sf2.writef(modulated.get_samples(), def);
    return 0;
}
