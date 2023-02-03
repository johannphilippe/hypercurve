#include <iostream>
#include<memory.h>
#include"../src/hypercurve.h"
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
    png p(2048, 1024, {255,255,255,0}, {244, 101, 36, 255});
    p.draw_curve(c.get_samples(), c.get_definition(), false, 0 , waveform);
    //p.draw_grid(10, 10, color{255, 255, 255, 100});
    std::string concat("/home/johann/Documents/" + name);
    p.write_as_png(concat);
}

void write_doc_png(curve &c, std::string name, int width = 1024, int height = 256)
{
    png p(width, height);
    p.draw_curve(c.get_samples(), c.get_definition(),true,0, false);
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
                    segment(fraction(2, segs), 0.8, share(mouth_curve())),
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
    spl.norm();
    write_doc_png(spl, "cubic_spline.png");
    curve catmul(2048, 0, {segment(1, 1, share(catmull_rom_spline_curve(point(-1,-0.5), point(2, 3.5))))});
    write_doc_png(catmul, "catmul_rom.png");
    curve polynomial(2048, 0, { segment(1, 1, share(polynomial_curve( {1.34, -1, -0.5, 0.1}) ))});
    polynomial.norm();
    write_doc_png(polynomial, "polynomial.png");
    curve typed(2048, 0, {segment(1, 1, share(typed_curve(-5)))});
    write_doc_png(typed, "typed.png");
    curve mouse(2048, 0, {segment(1, 1, share(mouth_curve()))});
    write_doc_png(mouse, "mouse.png");
    curve bicorn(2048, 0, {segment(1, 1, share(bicorn_curve(true)))});
    write_doc_png(bicorn, "bicorn.png");
    curve lagrange(2048, 0, {segment(1, 1, share(lagrange_polynomial_curve({ control_point(0.2, 0.8),  control_point(0.4, 0.1) })))});
    lagrange.scale(0, 1);
    write_doc_png(lagrange, "lagrange.png");

    curve log(2048, 0, {segment(1, 1, share(logarithmic_curve()))});
    write_doc_png(log, "logarithmic.png");
    curve exp(2048, 0, {segment(1, 1, share(exponential_curve()))});
    write_doc_png(exp, "exponential.png");
}


void unit_tests()
{
    int def = 16384;


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
                        point(0.3, 0.2),
                        point(0.6, 0.9),
                        point(0.8, 1)
                     }))
                 })
             });
    t.time_since_last("Cubic spline");
    c6.ascii_display("Cubic spline", "y = cspline(X)", '*');
    c6.scale(0.0, 1.0);
    write_as_png(c6,  false, "cubic_spline.png");

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

    modulated.scale(0, 1);
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
                 segment(fraction(1, hdiv), 0.5, share(cissoid_curve(0.51))),
                 segment(fraction(1, hdiv), 0.3, share(power_curve(9))),
                 segment(fraction(6, hdiv), 0, share(cubic_bezier_curve(point(0.2, 0.9), point(0.8, 0.8))))
              });

    c13.ascii_display("awesome hybrid curve", "A combination of gaussian, cissoid, power of 9, and cubic bezier curve", '*');


    curve c14(def, 0,
    {
                segment(fraction(1,2), 1, share(diocles_curve(1))),
                segment(fraction(1,2), 0, share(diocles_curve(1)))
              });

    c14.ascii_display("Diocles", "diocles with 1", 'y');

    curve cheb(def, 0,
    {
                  segment(1, 1, share( chebyshev_modulator<amplitude_fixed>(1, 20)))
               });

    cheb.ascii_display("chebyshev", "cheby", 'o');

    curve cat(def, 0,
    {
                 segment(1, 1, share(catenary_curve(0.1)))
              });
    cat.ascii_display("catenary " , "experimental...", '@');

    //write_as_png(cat, false, "catenary_a1.png");
    curve cat2(def, 0,
    {
                 segment(1, 1, share(catenary_curve(100000)))
              });

    cat2.ascii_display("catenary2 " , "experimental...", '@');
    //write_as_png(cat2, false, "catenary_a100");


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
                       segment(fraction(1,2), 1, vinvert( share(tightrope_walker_curve(1.01,1)))),
                       segment(fraction(1,2), 0, vinvert( share(tightrope_walker_curve(1.01,1)))),
                    });

    tightrope2.scale(-1, 1);
    tightrope2.ascii_display("tightrope walker2q", "trw", '*');
    write_as_png(tightrope2, true, "tightrope2.png");
    check_equality(tightrope, tightrope2);

    curve polynomial(def, 0, {
                         segment(1, 1, share(polynomial_curve( {1.34, -1, -0.5, 0.1}) ))
                     });

    polynomial.ascii_display("polynomial", "poly", '*');
    polynomial.scale(0, 1);
    write_as_png(polynomial, false, "poly.png");

    curve mouse(def, 0, {
                    segment(fraction(1,2), 1, share(mouth_curve())),
                    segment(fraction(1,2), 0, share(mouth_curve()))
                });
    mouse.ascii_display("mouse", "kiss", '*');
    write_as_png(mouse, false, "mouse.png");

    curve bicorn(def, 0, {
                     segment(fraction(1,2), 1, share(bicorn_curve(false))),
                     segment(fraction(1,2), 0, share(bicorn_curve(false))),
                 });
    bicorn.ascii_display("bicorn", "cocked hat", '*');
    bicorn.scale(0, 1);
    write_as_png(bicorn, false, "bicorn.png");

    // Curves pictures for Doc
    curve rescale_test(def, 0, {
                           segment(fraction(1, 3), 1, share(bicorn_curve(false))),
                           segment(fraction(1, 3), 0, share(mouth_curve()))
                       });
    rescale_test.ascii_display("Rescaled", "rescaled 2/3", '*');



    curve log(def, 0, {
                 segment(fraction(1, 2), 1, share(logarithmic_curve())),
                 segment(fraction(1, 2), 0, share(logarithmic_curve())),
              });

    log.ascii_display("logarithmic scaled", "not rescaled", '*');
    write_as_png(log, false, "log.png");

    // exponential as mirrored log
    curve xp(def, 0, {
                 segment(fraction(1, 2), 1, mirror(share(logarithmic_curve()))),
                 segment(fraction(1, 2), 0, mirror(share(logarithmic_curve()))),
              });

    xp.ascii_display("exponential scaled", "not rescaled", '*');
    write_as_png(xp, false, "exp.png");

    curve cmp(def, 0, {
                 segment(fraction(1, 2), 1, share(power_curve(2))),
                 segment(fraction(1, 2), 0, share(power_curve(2))),

              });
    cmp.ascii_display("cmp power", "not rescaled", '*');

    check_equality(xp, cmp);

    curve logarithm(def, 0, {
                       segment(1, 1, share(logarithmic_curve()))
                    });
    write_as_png(logarithm, false, "LOGARITHM.png");



}

void random_generator_test()
{
    // Random curves
    const int def = 16384;

    for(size_t i = 0; i < 100; ++i){
        std::pair<curve, std::string> rand1 = random_curve_composer(4, 0, 1, def, true, false);
        rand1.first.ascii_display("random_curve" + std::to_string(i), rand1.second, '*');

        std::string rname = "rand_" + std::to_string(i) + "__" + rand1.second + ".png";
        write_as_png(rand1.first, false, rname);
    }

}

void bark_bands_test()
{
    size_t def = 16384;
    size_t hdiv = 9;
    curve c13(def, 0, {
                 segment(fraction(1, hdiv), 1, share(gauss_curve(5, 1))),
                 segment(fraction(1, hdiv), 0.5, share(cissoid_curve(0.51))),
                 segment(fraction(1, hdiv), 0.3, share(power_curve(9))),
                 segment(fraction(6, hdiv), 0, share(cubic_bezier_curve(point(0.2, 0.9), point(0.8, 0.8))))
              });

    c13.ascii_display("Not BARK", "NOOOOOO", '*');
    write_as_png(c13, false, "not_bark.png");
    scale_to_bark(c13.get_samples(), def);
    c13.ascii_display("BARK", "BAAAAARK", '*');
    c13.scale(0, 1);
    write_as_png(c13, false, "bark.png");
}

void dummy_test()
{

    size_t def = 16384;


    // Tests on png

    curve xplode(def, 0, {
                     segment(1, 1, share(lagrange_polynomial_curve({point(0.01, 0.9), point(0.02, 0.01), point(0.021, 0.99)})))
                 });
    xplode.ascii_display("xplode", "lagrange polynomial", '*');
    write_as_png(xplode, false, "xplode.png");


    curve wav(def, 0, {
                  segment(fraction(1,3), 1, share(lagrange_polynomial_curve({point(0.01, 0.9), point(0.02, 0.01), point(0.021, 0.99)}))),
                  segment(fraction(1,3), -1, share(lagrange_polynomial_curve({point(0.01, 0.9), point(0.02, 0.01), point(0.021, 0.99)}))),
                  segment(fraction(1,3), 0, share(lagrange_polynomial_curve({point(0.01, 0.9), point(0.02, 0.01), point(0.021, 0.99)})))
              });
    wav.ascii_display("waveform", "lagrange wav", '*');
    write_as_png(wav, false, "waveform.png");
    /*
    curve cnorm(def, 0, {
                 segment(1, 1, (share(cubic_bezier_curve(point(0.1, 0.95), point(0.2, -1.01)))))
             });
    cnorm.ascii_display("bezier curve", "bez", '*');
    curve cinv(def, 0, {
                 segment(1, 1, hinvert(share(cubic_bezier_curve(point(0.1, 0.95), point(0.2, -1.01)))))
             });
    cinv.ascii_display("hinverted bezier curve", "bez", '*');
    curve cdinv(def, 0, {
                 segment(1, 1, vinvert(hinvert(share(cubic_bezier_curve(point(0.1, 0.95), point(0.2, -1.01))))))
             });
    cdinv.ascii_display("hinverted and vinverted bezier curve", "bez", '*');
    */

    /*
    curve c2norm(def, 0, {
                     segment(1, 1, (share(tightrope_walker_curve(1.01, 1))))
             });
    c2norm.ascii_display("bezier curve", "bez", '*');
    write_as_png(c2norm, false, "noinv.png");
    curve c2inv(def, 0, {
                 segment(1, 1, hinvert(share(tightrope_walker_curve(1.01, 1))))
             });
    c2inv.ascii_display("hinverted bezier curve", "bez", '*');
    write_as_png(c2inv, false, "hinv.png");

    curve c2vdinv(def, 0, {
                     segment(1, 1, vinvert((share(tightrope_walker_curve(1.01, 1)))))
             });
    c2vdinv.ascii_display(" vinverted bezier curve", "bez", '*');
    c2vdinv.norm();
    write_as_png(c2vdinv, false, "vinv.png");

    curve c2dinv(def, 0, {
                     segment(1, 1, vinvert(hinvert(share(tightrope_walker_curve(1.01, 1)))))
             });
    c2dinv.ascii_display("hinverted and vinverted bezier curve", "bez", '*');
    c2dinv.norm();
    write_as_png(c2dinv, false, "vhinv.png");
    */

    /*
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
    */
    /*
    curve crv(16384, 0, {segment(1., 1., share(
                       lagrange_polynomial_curve(
                       { control_point(0.2, 0.9),
                         control_point(0.4, 0.01),
                         control_point(0.7, 0.9)
                       }))
                       )});
    crv.ascii_display("lagrange curve", "interpolation", '*');
    crv.scale(0, 1);
    write_as_png(crv, false, "lagrange.png");

    curve crvspl(16384, 0, {segment(1., 1., share(
                       cubic_spline_curve(
                       { control_point(0.2, 0.9),
                         control_point(0.4, 0.01),
                         control_point(0.7, 0.9)
                       }))
                       )});
    crvspl.ascii_display("spline curve", "interpolation", '*');
    crvspl.scale(0, 1);
    write_as_png(crvspl, false, "spline.png");
    */

}


void inversion_test()
{
    const double start = 0.3;
    const double end = 0.8;
    curve dio(4096, start, {
                 segment(0.5, end, share(tightrope_walker_curve(1.1, 0.1))),
                 segment(0.5, 0, share(diocles_curve(2)))
              });
    dio.ascii_display("Diocles", "dio", '*');
    write_as_png(dio, false, "dio.png");

    curve inv(4096, start, {
                 segment(0.5, end, vinvert(share(tightrope_walker_curve(1.1, 0.1)))),
                 segment(0.5, 0, vinvert(share(diocles_curve(2))))
              });

    inv.ascii_display("Inverted", "Inv", '*');
    inv.scale(0, 1);
    write_as_png(inv, false, "inv.png");


    AsciiPlotter p("comp", 80, 30);
    std::vector<double> dio_vec(dio.get_definition());
    ::memcpy(dio_vec.data(), dio.get_samples(), sizeof(double) * dio.get_definition());
    std::vector<double> inv_vec(dio.get_definition());
    ::memcpy(inv_vec.data(), inv.get_samples(), sizeof(double) * inv.get_definition());


    curve linear_axis(4096, start, {
                         segment(0.5, end, share(linear_curve())),
                         segment(0.5, 0, share(linear_curve()))
                      });
    std::vector<double> linear_vec(linear_axis.get_definition());
    ::memcpy(linear_vec.data(), linear_axis.get_samples(), sizeof(double) * linear_axis.get_definition());

    p.addPlot(dio_vec, "diocles", '\\');
    p.addPlot(inv_vec, "inverted", '/');


    p.addPlot(linear_vec, "linear", '|');
    p.legend();
    p.show();
}


void icsc()
{
    double def = 8192;
    double div = 64;
    double atq = 3;
    double dec = 24;
    double sus = 0.35;
    double rel = div - (atq + dec);

    curve tightrope(def, 0, {
                       segment(fraction(atq, div), 1, share(tightrope_walker_curve(1.105, 0.125))),
                       segment(fraction(dec, div), sus, share(tightrope_walker_curve(0.95, 0.25))),
                       segment(fraction(rel, div), 0, share(tightrope_walker_curve(0.5, 0.15)))
                    });
    write_as_png(tightrope, false, "icsc_tightrope.png");
    curve kiss(def, 0, {
                       segment(fraction(atq, div), 1, share(kiss_curve())),
                       segment(fraction(dec, div), sus, share(kiss_curve())),
                       segment(fraction(rel, div), 0, share(kiss_curve()))
                    });
    write_as_png(kiss, false, "icsc_kiss.png");

    curve cate(def, 0, {
                       segment(fraction(atq, div), 1, share(catenary_curve(1.75))),
                       segment(fraction(dec, div), sus, share(catenary_curve(0.95))),
                       segment(fraction(rel, div), 0, share(catenary_curve(0.25)))
                    });
    write_as_png(cate, false, "icsc_catenary.png");
    curve tox(def, 0, {
                       segment(fraction(atq, div), 1, share(toxoid_curve(0.05))),
                       segment(fraction(dec, div), sus, share(toxoid_curve(5.95))),
                       segment(fraction(rel, div), 0, share(toxoid_curve(0.05)))
                    });
    write_as_png(tox, false, "icsc_toxoid.png");


    curve simple(def, 0, {
                     segment(0.1, 1, share(diocles_curve(0.5001))),
                     segment(0.9, 0, share(diocles_curve(0.55)))
                 });
    write_as_png(simple, false, "icsc_diocles.png");


    curve complex(def, 0, {
                      segment(0.05, 1, share(gauss_curve(1, 0.1))),
                      segment(0.4, 0.5, share(catmull_rom_spline_curve(point(-1, 3), point(3, -2)))),
                      segment(0.2, 0.3, mirror(share(bicorn_curve(true)))),
                      segment(0.35, 0, share(cubic_bezier_curve(point(0.1, 0.1), point(0.9, 0.65))))
                      //segment(0.3, 0.0, share(linear_curve()))
                  });

    //complex.norm();
    complex.find_extremeness();
    complex.ascii_display("complex", "heeey", '*');
    write_as_png(complex, false, "icsc_complex.png");
    /*
    for(size_t i = 0; i < 100; ++i) {
        std::pair<curve, std::string> rnd = random_curve_composer(6, 0, 1, 16384, true, false, false );
        std::string name = "purpleicsc_" + rnd.second + ".png";
        write_as_png(rnd.first, false, name);

    }
    */
}

int main()
{

    //icsc();
    //dummy_test();
    //unit_tests();
    //random_generator_test();
    //inversion_test();
    generate_curve_pictures();

    /*
    auto crv = curve(16384, 0, {
                        segment(0.5, 1, share(diocles_curve(0.51))),
                        segment(0.5, 0, share(tightrope_walker_curve(1.1, 0.1)))
                     });
    write_as_png(crv, false, "dio_tightrope.png");

    auto crv2 = curve(16384, 0, {
                         segment(0.1, 1, share(tightrope_walker_curve(1.1, 0.1))),
                         segment(0.4, 0.2, share(gauss_curve(1, 1))),
                         segment(0.5, 0, share(kiss_curve()))
                      });
    write_as_png(crv2, false, "tightrope_gauss_kiss.png");
    */


    return 0;
}
