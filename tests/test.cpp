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
    SndfileHandle sf("test_hypercurve.wav", SFM_WRITE, fmt, 1, 48000);
    // composite curve
    curve c(def, 0
            , 	{
                share(segment(frac(1,4), 1.0, share(diocles_curve(1)))),
                share(segment(frac(1,4), 0.5, share(cubic_curve()))),
                share(segment(frac(1,4), 1, share(diocles_curve(1)))),
                share(segment(frac(1,4), 0, share(diocles_curve(1))))
                }
            );
    // One segment curve
    curve c2(def, 0, {share(segment(frac(1,1), 1.0, share(diocles_curve(1))))});
    double div = 4;

    // Bezier quadratic (one control point)
    curve c3(def, 0, {
                 share(bezier_segment(frac(1,div), 1.0, share(quadratic_bezier_curve({0.1, 0.9})))),
                 share(bezier_segment(frac(1,div), 0.5, share(quadratic_bezier_curve({0.66, 0.1})))),
                 share(bezier_segment(frac(1,div), 0.8, share(quadratic_bezier_curve({0.9, 0.9})))),
                 share(bezier_segment(frac(1, div), 0.1, share(quadratic_bezier_curve({0.5, 0}))))
             });

    // Bezier cubic (taking two control points)
    curve c4(def, 0, {
                 share(bezier_segment(frac(1,1), 1	, share(cubic_bezier_curve({0.2, 0.8}, {0.8, 0.2}))))
             });

    // homemade smooth
    curve c5(def, 0, {
                share(segment(1, 1, share(hypersmooth_curve())))
             });

    sf.writef(c4.get_samples(), def);

    return 0;
}
