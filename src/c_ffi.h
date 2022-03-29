#ifndef C_FFI_H
#define C_FFI_H
#include"curve_lib.h"
#include"core.h"

extern "C" {

    /*
     *
     * CORE
     *
    */

    segment c_segment(double frac, double y_dest, std::shared_ptr<curve_base> crv)
    {
        return segment(frac, y_dest, crv);
    }


    // Find a way to build a segment vector from Ruby/Lua
    std::vector<segment> c_segs()
    {

    }

    // And a way to return a Lua/Ruby array full of pretty curved samples
    curve *c_curve(size_t def, double y_start, segment segs[])
    {
        return new curve(def, y_start, std::vector<segment>(segs, ) );
    }

    /*
     *
     * CURVE LIBRARY
     *
    */

    std::shared_ptr<diocles_curve> c_diocles_curve(double a_)
    {
        return share(diocles_curve(a_));
    }

    std::shared_ptr<cubic_curve> c_cubic_curve()
    {
        return share(cubic_curve());
    }





}

#endif // C_FFI_H
