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
    p.draw_grid(10, 10, color{255, 255, 255, 100});
    std::string concat("/home/johann/Documents/" + name);
    p.write_as_png(concat);
}

void write_doc_png(curve &c, std::string name, int width = 1024, int height = 256)
{
    png p(width, height);
    p.draw_curve(c.get_samples(), c.get_definition(), true, false);
    p.draw_grid(10, 10, color{255, 255, 255, 100});
    std::string concat("/home/johann/Documents/GitHub/hypercurve/doc/png/" + name);
    std::cout << concat << std::endl;
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


void generate_curve_pictures()
{
    const int segs = 8;
    curve hybrid(2048, 0, {
                    segment(fraction(1, segs), 1, share(diocles_curve(1))),
                    segment(fraction(1, segs), 0.2, share(toxoid_curve(0.01))),
                    segment(fraction(2, segs), 0.8, share(mouse_curve())),
                    segment(fraction(4, segs), 0, share(gauss_curve(10, 0.5)))
                 });
    hybrid.ascii_display("hybrid_curve", "hybrid", '*');
    write_doc_png(hybrid,  "hybrid.png", 2048, 512);

    curve dioc(2048, 0, {segment(1, 1, share(diocles_curve(1)))});
    write_doc_png(dioc, "diocles.png");
    curve cubic(2048, 0, {segment(1, 1, share(cubic_curve()))});
    write_doc_png(cubic, "cubic.png");
    curve power(2048, 0, {segment(1, 1, share(power_curve(9)))});
    write_doc_png(power, "power9.png");
    curve hann(2048, 0, {segment(1, 1, share(hanning_curve()))});
    write_doc_png(hann, "hanning.png");
    curve hamm(2048, 0, {segment(1, 1, share(hamming_curve()))});
    write_doc_png(hamm, "hamming.png");
    curve black(2048, 0, {segment(1, 1, share(blackman_curve()))});
    write_doc_png(black, "blackman.png");
    curve gauss(2048, 0, {segment(1, 1, share(gauss_curve(10, 0.5)))});
    write_doc_png(gauss, "gaussian.png");
    curve tox(2048, 0, {segment(1, 1, share(toxoid_curve(10)))});
    write_doc_png(tox, "toxoid.png");
    curve cat(2048, 0, {segment(1, 1, share(catenary_curve(0.1)))});
    write_doc_png(cat, "catenary.png");
    curve tight(2048, 0, {segment(1, 1, share(tightrope_walker_curve(1.1, 0.1)))});
    write_doc_png(tight, "tightrope.png");
    curve quadbez(2048, 0, {segment(1, 1, share(quadratic_bezier_curve(point(0.2, 0.1))))});
    write_doc_png(quadbez, "quadratic_bezier.png");
    curve cubbez(2048, 0, {segment(1, 1, share(cubic_bezier_curve(point(0.1, 0.1), point(0.5, 0.8))))});
    write_doc_png(cubbez, "cubic_bezier.png");
    curve spl(2048, 0, {segment(1, 1, share(cubic_spline_curve( { point(0.2, 0.7), point(0.3, 0.2), point(0.5, 0.8)} )))});
    write_doc_png(spl, "cubic_spline.png");
    curve catmul(2048, 0, {segment(1, 1, share(catmull_rom_spline_curve(point(-1,-0.5), point(2, 3.5))))});
    write_doc_png(catmul, "catmul_rom.png");
    curve polynomial(2048, 0, { segment(1, 1, share(polynomial_curve( {1.34, -1, -0.5, 0.1}) ))});
    polynomial.normalize_y(0, 1);
    write_doc_png(polynomial, "polynomial.png");
    curve typed(2048, 0, {segment(1, 1, share(typed_curve(-5)))});
    write_doc_png(typed, "typed.png");
    curve mouse(2048, 0, {segment(1, 1, share(mouse_curve()))});
    write_doc_png(mouse, "mouse.png");
    curve bicorn(2048, 0, {segment(1, 1, share(bicorn_curve(true)))});
    write_doc_png(bicorn, "bicorn.png");
    curve lagrange(2048, 0, {segment(1, 1, share(lagrange_polynomial_curve({ control_point(0.2, 0.8),  control_point(0.4, 0.1) })))});
    lagrange.normalize_y(0, 1);
    write_doc_png(lagrange, "lagrange.png");
}


void unit_tests()
{
    int def = 16384;

    const int fmt = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
    SndfileHandle sf("test_hypercurve_HYBRID.wav", SFM_WRITE, fmt, 1, 48000);
    SndfileHandle sf2("test_hypercurveGauss2.wav", SFM_WRITE, fmt, 1, 48000);

    timer t;

    // Simple composite curve
    curve c(def, 0
            , 	{
                segment(fraction(1,4), 1.0, share(diocles_curve(1))),
                segment(fraction(1,4), 0.5, share(cubic_curve())),
                segment(fraction(1,4), 1, share(diocles_curve(1))),
                segment(fraction(1,4), 0, share(diocles_curve(1)))
                }
            );

    t.time_since_last("composite curve");

    c.ascii_display("Composite diocles and cubic curve", "y = composite(x)", '*');

    t.reset();

    // One segment curve
    curve c2(def, 0, {segment(fraction(1,1), 1.0, share(diocles_curve(1)))});

    t.time_since_last("one segment curve");

    c2.ascii_display("One segment diocles curve", "y = oneseg(x)", '*');

    const double div = 4;
    t.reset();
    // Bezier quadratic (one control point)
    curve c3(def, 0, {
                 segment(fraction(1,div), 1.0, share(quadratic_bezier_curve({0.1, 0.9}))),
                 segment(fraction(1,div), 0.5, share(quadratic_bezier_curve({0.66, 0.1}))),
                 segment(fraction(1,div), 0.8, share(quadratic_bezier_curve({0.5, 0.9}))),
                 segment(fraction(1, div), 0.1, share(quadratic_bezier_curve({0.5, 0})))
             });

    t.time_since_last("Bezier quadratic 4 segments");
    write_as_png(c3, false, "bezier_quadratic.png");

    c3.ascii_display("Bezier quadratic", "y = quadratic_bezier(x)", '*');

    t.reset();
    // Bezier cubic (taking two control points)
    curve c4(def, 0, {
                 segment(fraction(1,1), 1	, share(cubic_bezier_curve({0.2, 0.8}, {0.8, 0.2})))
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
                 segment(fraction(1,2), 1, {
                     share(catmull_rom_spline_curve(
                        point(-2, -0.5),
                        point(2, 0.2)))
                 }),
                 segment(fraction(1,2), 0, {
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
                  segment(fraction(1,3), 1, share(cubic_curve())),
                  segment(fraction(1,3), -1, share(cubic_curve())),
                  segment(fraction(1,3), 0, share(cubic_curve())),

               });
    t.time_since_last("Waveform");

    w.ascii_display("waveform", "w=cubic(x)", '.');


    curve to_modulate(16384, 0, {
                         segment(fraction(1,2), 1, share(diocles_curve(1))),
                         segment(fraction(1,2), 0, share(diocles_curve(1))),

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
                 segment(fraction(1,2), 1, share(hanning_curve())),
                 segment(fraction(1,2), 0, share(hanning_curve()))
             });
    c8.ascii_display("Hanning", "hanning(x)", '.');


    curve c9(def, 0.0, {
                segment(1, 1, share(gauss_curve(10, 0.5	)))
             });
    c9.ascii_display("gauss", "A = 10, c = 0.5", '*');
    curve c10(def, 0.0, {
                segment(fraction(1,2), 1, share(gauss_curve(1, 0.5))),
                segment(fraction(1,2), 0, share(gauss_curve(1, 0.5)))
             });
    c10.ascii_display("gauss", "A = 1, c = 0.5", '*');
    write_as_png(c10, false, "gauss.png");

    check_equality(c9, c10);

    curve c11(def, 0, {
                 segment(fraction(1,2), 1, share(typed_curve(10))),
                 segment(fraction(1,2), 0, share(typed_curve(10)))
              });
    c11.ascii_display("typed curve" , "type = 10", '*');
    write_as_png(c11, false, "typedpos.png");
    curve c12(def, 0, {
                 segment(fraction(1,2), 1, share(typed_curve(-10))),
                 segment(fraction(1,2), 0, share(typed_curve(-10)))
              });
    c12.ascii_display("typed curve" , "type = -10", '*');
    write_as_png(c12, false, "typedneg.png");

    double hdiv = 9;
    curve c13(def, 0, {
                 segment(fraction(1, hdiv), 1, share(gauss_curve(5, 1))),
                 segment(fraction(1, hdiv), 0.8, share(cissoid_curve(2))),
                 segment(fraction(1, hdiv), 0.1, share(power_curve(9))),
                 segment(fraction(6, hdiv), 0, share(cubic_bezier_curve(point(0.2, 0.9), point(0.8, 0.8))))
                 //segment(fraction(6, hdiv), 0, share(cissoid_curve(1)))
              });

    c13.ascii_display("awesome hybrid curve", "A combination of gaussian, cissoid, power of 9, and cubic bezier curve", '*');


    curve c14(def, 0,
    {
                segment(fraction(1,2), 1, share(diocles_curve(1))),
                segment(fraction(1,2), 0, share(diocles_curve(1)))
              });

    curve cheb(def, 0,
    {
                  segment(1, 1, share( chebyshev_modulator<amplitude_fixed>(1, 20)))
               });



    curve cat(def, 0,
    {
                 segment(1, 1, share(catenary_curve(0.1)))
              });

    write_as_png(cat, false, "catenary_a1.png");
    curve cat2(def, 0,
    {
                 segment(1, 1, share(catenary_curve(100000)))
              });

    write_as_png(cat2, false, "catenary_a100");


    curve tox(def, 0, {
                 segment(fraction(1, 2), 1, share(toxoid_curve(10))),
                 segment(fraction(1, 2), 0, share(toxoid_curve(0.5)))
              });
    tox.ascii_display("toxoid", "tox", '*');
    write_as_png(tox, false, "toxoid.png");

    curve ctst(def, 0, {
                       segment(fraction(1,10), 1, share(power_curve(50))),
                       segment(fraction(9, 10), 0, share(power_curve(9)))
                    });
    ctst.ascii_display("ctst", "power curve", '*');
    write_as_png(ctst, false, "power.png");

    curve tightrope(def, 0, {
                       segment(fraction(2,10), 1, share(tightrope_walker_curve(1.1,0.1))),
                       segment(fraction(8, 10), 0, share(tightrope_walker_curve(1.1, 1.0)))
                    });
    tightrope.ascii_display("tightrope walker", "trw", '*');
    write_as_png(tightrope, false, "tightrope.png");

    curve tightrope2(def, 0, {
                       segment(fraction(1,2), 1, invert( share(tightrope_walker_curve(1.01,1)))),
                       segment(fraction(1,2), 0, invert( share(tightrope_walker_curve(1.01,1)))),
                    });

    tightrope2.normalize_y(-1, 1);
    tightrope2.ascii_display("tightrope walker2q", "trw", '*');
    write_as_png(tightrope2, true, "tightrope2.png");
    check_equality(tightrope, tightrope2);

    sf.writef(c13.get_samples(), def);
    sf2.writef(c10.get_samples(), def);


    curve polynomial(def, 0, {
                         segment(1, 1, share(polynomial_curve( {1.34, -1, -0.5, 0.1}) ))
                     });

    polynomial.ascii_display("polynomial", "poly", '*');
    polynomial.normalize_y(0, 1);
    write_as_png(polynomial, false, "poly.png");

    curve mouse(def, 0, {
                    segment(fraction(1,2), 1, share(mouse_curve())),
                    segment(fraction(1,2), 0, share(mouse_curve()))
                });
    mouse.ascii_display("mouse", "kiss", '*');
    write_as_png(mouse, false, "mouse.png");

    curve bicorn(def, 0, {
                     segment(fraction(1,2), 1, share(bicorn_curve(false))),
                     segment(fraction(1,2), 0, share(bicorn_curve(false))),
                 });
    bicorn.ascii_display("bicorn", "cocked hat", '*');
    bicorn.normalize_y(0, 1);
    write_as_png(bicorn, false, "bicorn.png");

    // Curves pictures for Doc
    curve rescale_test(def, 0, {
                           segment(fraction(1, 3), 1, share(bicorn_curve(false))),
                           segment(fraction(1, 3), 0, share(mouse_curve()))
                       });
    rescale_test.ascii_display("Rescaled", "rescaled 2/3", '*');


}

void random_generator_test()
{
    // Random curves
    const int def = 16384;

    for(size_t i = 0; i < 100; ++i){
        std::pair<curve, std::string> rand1 = random_curve_composer(4, -1, 1, def, true, true);
        rand1.first.ascii_display("random_curve" + std::to_string(i), rand1.second, '*');

        std::string rname = "rand_" + std::to_string(i) + "__" + rand1.second + ".png";
        write_as_png(rand1.first, true, rname);
    }

}

void dummy_test()
{

    curve crv(16384, 0, {segment(0.5, 1, share(hamming_curve())), segment(0.5, 0, share(hamming_curve()))});
    write_as_png(crv, false, "hamming.png");
    for(size_t i = 0; i < 1024; ++i)
    {
        std::cout << "hanning at "  << i << " is : " << hamming(i, 1024) << std::endl;
    }
    curve crvb(16384, 0, {segment(fraction(1,2), 1, share(blackman_curve())), segment(fraction(1,2), 0, share(blackman_curve()))});
    write_as_png(crvb, false, "blackman.png");
    for(size_t i = 0; i < 1024; ++i)
    {

        std::cout << "blackman at "  << i << " is : " << blackman(i, 1024) << std::endl;
    }
    /*
    curve crv(16384, 0, {segment(1., 1., share(
                       lagrange_polynomial_curve(
                       { control_point(0.2, 0.9),
                         control_point(0.4, 0.01),
                         control_point(0.7, 0.9)
                       }))
                       )});
    crv.ascii_display("lagrange curve", "interpolation", '*');
    crv.normalize_y(0, 1);
    write_as_png(crv, false, "lagrange.png");

    curve crvspl(16384, 0, {segment(1., 1., share(
                       cubic_spline_curve(
                       { control_point(0.2, 0.9),
                         control_point(0.4, 0.01),
                         control_point(0.7, 0.9)
                       }))
                       )});
    crvspl.ascii_display("spline curve", "interpolation", '*');
    crvspl.normalize_y(0, 1);
    write_as_png(crvspl, false, "spline.png");
    */

}

int main()
{

    //dummy_test();
    //unit_tests();
    random_generator_test();
    //generate_curve_pictures();

    return 0;
}
