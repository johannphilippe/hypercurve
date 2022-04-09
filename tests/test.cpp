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
void write_as_png(curve &c, bool waveform = true, std::string name = "h.png")
{
    png p(2048, 1024);
    p.draw_curve(c.get_samples(), c.get_definition(), true, waveform);
    p.draw_grid(10, 10, color{{255, 255, 255, 100}});
    std::string concat("/home/johann/Documents/" + name);
    p.write_as_png(concat);
}

void check_equality(curve &c1, curve &c2)
{
    for(size_t i = 0; i < c1.get_definition(); i++)
    {
        if(c1.get_sample_at(i) != c2.get_sample_at(i))
        {
            std::cout << "ARE NOT IDENTICAL" << std::endl;
            return;
        }
    }
    std::cout << "IDENTICAL "  << std::endl;

}

int main()
{
    int def = 16384;

    const int fmt = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
    SndfileHandle sf("test_hypercurve_HYBRID.wav", SFM_WRITE, fmt, 1, 48000);
    SndfileHandle sf2("test_hypercurveGauss2.wav", SFM_WRITE, fmt, 1, 48000);

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
                segment(1, 1, share(power_curve(9)))
             });

    t.time_since_last("Homemade power of 9");
    c5.ascii_display("Power curve", "y = x^9", '*');

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

    cubic_interpolator itp(0, { point(0.5, 0.2), point(0.5, 0)});
    curve mod(16384, 0, {
                  segment(1, 1, share(chebyshev_modulator<amplitude_interpolated>(share(itp), 15)))
            });

    mod.ascii_display("modulator" ,"mod(x)", '*');
    curve modulated = to_modulate + mod;

    to_modulate.ascii_display("to_modulate", "tomod", '*');

    modulated.normalize_y(0, 1);
    modulated.ascii_display("Modulated", "mod * catmullrom", '*');


    curve c8(def, 0.0, {
                 segment(frac(1,2), 1, share(hanning_curve())),
                 segment(frac(1,2), 0, share(hanning_curve()))
             });
    c8.ascii_display("Hanning", "hanning(x)", '.');


    curve c9(def, 0.0, {
                segment(1, 1, share(gauss_curve(10, 0.5	)))
             });
    c9.ascii_display("gauss", "A = 10, c = 0.5", '*');
    curve c10(def, 0.0, {
                segment(1, 1, share(gauss_curve(1, 0.5)))
             });
    c10.ascii_display("gauss", "A = 1, c = 0.5", '*');

    check_equality(c9, c10);

    curve c11(def, 0, {
                 segment(frac(1,2), 1, share(typed_curve(10))),
                 segment(frac(1,2), 0, share(typed_curve(10)))
              });
    c11.ascii_display("typed curve" , "type = 10", '*');
    curve c12(def, 0, {
                 segment(frac(1,2), 1, share(typed_curve(-10))),
                 segment(frac(1,2), 0, share(typed_curve(-10)))
              });
    c12.ascii_display("typed curve" , "type = -10", '*');

    double hdiv = 9;
    curve c13(def, 0, {
                 segment(frac(1, hdiv), 1, share(gauss_curve(5, 1))),
                 segment(frac(1, hdiv), 0.8, share(cissoid_curve(2))),
                 segment(frac(1, hdiv), 0.1, share(power_curve(9))),
                 segment(frac(6, hdiv), 0, share(cubic_bezier_curve(point(0.2, 0.9), point(0.8, 0.8))))
                 //segment(frac(6, hdiv), 0, share(cissoid_curve(1)))
              });

    c13.ascii_display("awesome hybrid curve", "A combination of gaussian, cissoid, power of 9, and cubic bezier curve", '*');


    curve c14(def, 0,
    {
                segment(frac(1,2), 1, share(diocles_curve(1))),
                segment(frac(1,2), 0, share(diocles_curve(1)))
              });

    curve cheb(def, 0,
    {
                  segment(1, 1, share( chebyshev_modulator<amplitude_fixed>(1, 20)))
               });



    curve cat(def, 0,
    {
                 segment(1, 1, share(catenary_curve(0.1)))
              });

    write_as_png(cat, false, "catenary_a1");
    curve cat2(def, 0,
    {
                 segment(1, 1, share(catenary_curve(100000)))
              });

    write_as_png(cat2, false, "catenary_a100");


    curve tox(def, 0, {
                 segment(frac(1, 2), 1, share(toxoid_curve(10))),
                 segment(frac(1, 2), 0, share(toxoid_curve(0)))
              });
    write_as_png(tox, false, "toxoid");

    sf.writef(c13.get_samples(), def);
    sf2.writef(c10.get_samples(), def);
    return 0;
}
