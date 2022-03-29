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
    curve c(def, 0
            , 	{
                segment(frac(1,4), 1.0, share(diocles_curve(1))),
                segment(frac(1,4), 0.5, share(cubic_curve())),
                segment(frac(1,4), 1, share(diocles_curve(1))),
                segment(frac(1,4), 0, share(diocles_curve(1)))
                }
            );

    curve c2(def, 0, {segment(frac(1,1), 1.0, share(diocles_curve(1)))});


    sf.writef(c.get_samples(), def);

    return 0;
}
