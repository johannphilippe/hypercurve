#include <iostream>
#include<sndfile.hh>
#include<memory.h>
#include"../src/core.h"
#include"../src/curve_lib.h"
#include"sndfile.hh"

using namespace std;

/*
 *  Call curve with full size
 *
 * 		// definition 	// starty
 * curve(16384, 		0.1,
 * {		// size		//toy	// mode (default to linear)
 * 	segment(4/8,	, 0.5	 	cubic_curve()),
 * });
*/

using namespace hypercurve;

int main()
{
    int def = 16384;

    const int fmt = SF_FORMAT_WAV | SF_FORMAT_PCM_32;
    SndfileHandle sf("test_hypercurveChebyshev1.wav", SFM_WRITE, fmt, 1, 48000);
    SndfileHandle sf2("test_hypercurveChebyshev2.wav", SFM_WRITE, fmt, 1, 48000);

    // Simple composite curve
    curve c(def, 0
            , 	{
                segment(frac(1,4), 1.0, share(diocles_curve(1))),
                segment(frac(1,4), 0.5, share(cubic_curve())),
                segment(frac(1,4), 1, share(diocles_curve(1))),
                segment(frac(1,4), 0, share(diocles_curve(1)))
                }
            );
    c.ascii_display("Composite diocles and cubic curve", "y = composite(x)", '*');

    // One segment curve
    curve c2(def, 0, {segment(frac(1,1), 1.0, share(diocles_curve(1)))});
    c2.ascii_display("One segment diocles curve", "y = oneseg(x)", '*');
    double div = 4;

    // Bezier quadratic (one control point)
    curve c3(def, 0, {
                 segment(frac(1,div), 1.0, share(quadratic_bezier_curve({0.1, 0.9}))),
                 segment(frac(1,div), 0.5, share(quadratic_bezier_curve({0.66, 0.1}))),
                 segment(frac(1,div), 0.8, share(quadratic_bezier_curve({0.9, 0.9}))),
                 segment(frac(1, div), 0.1, share(quadratic_bezier_curve({0.5, 0})))
             });
    c3.ascii_display("Bezier quadratic", "y = quadratic_bezier(x)", '*');

    // Bezier cubic (taking two control points)
    curve c4(def, 0, {
                 segment(frac(1,1), 1	, share(cubic_bezier_curve({0.2, 0.8}, {0.8, 0.2})))
             });
    c4.ascii_display("Bezier Cubic", "y = cubic_bezier(x)", '*');

    // homemade smooth
    curve c5(def, 0, {
                segment(1, 1, share(exponent_curve(9)))
             });
    c5.ascii_display("Smooth homemade", "y = hsmooth(x)", '*');

    // Cubic spline
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
    c6.ascii_display("Cubic spline", "y = cspline(X)", '*');

    // CatmullRom
    curve c7(def, 0, {
                 segment(frac(1,2), 1, {
                     share(catmull_rom_spline_curve(0.5,
                        point(-2, -0.5),
                        point(2, 0.2)))
                 }),
                 segment(frac(1,2), 0, {
                     share(catmull_rom_spline_curve(0.5,
                     point(-1, 3),
                     point(3, -5)
                     ))
                 })
             });

    c7.ascii_display("Composite", "y = composite(X)", '*');

    waveform w(def, 0, {
                  segment(frac(1,3), 1, share(cubic_curve())),
                  segment(frac(1,3), -1, share(cubic_curve())),
                  segment(frac(1,3), 0, share(cubic_curve())),

               });

    w.ascii_display("waveform", "w=cubic(x)", '.');

    const int n = 100;
    curve c8(def, 0, {
                segment(frac(1,1), 1, share(chebyshev_curve<1>(n)))
             });

    c8.ascii_display("Chebyshev type 1 n=10", "y=cheb(x)", '*');

    curve c9(def, 0, {
                segment(frac(1,1), 1, share(chebyshev_curve<2>(n)))
             });

    c9.ascii_display("Chebyshev type 2 n=10", "y=cheb(x)", '*');
    sf.writef(c8.get_samples(), def);
    sf2.writef(c9.get_samples(), def);

    return 0;
}
