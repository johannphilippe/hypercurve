# Hypercurve documentation

Hypercurve is a library of 2D curves designed to process audio envelopes, applied to any audio parameter. 
It is available in several frontends : C++, Lua, and Csound. 


## Hypercurve syntax

Here is a simple example of syntax with possible use cases : 

C++ : 
```c++
    auto crv = hypercurve::hypercurve( 2048, 0, {
        hypercurve::segment(0.5, 1, cubic_curve()),
        hypercurve::segment(0.5, 0, diocles_curve(1))
    });

    double sample = crv.get_sample_at(1024);
    double *samples = crv.get_samples();

    hypercurve::png png;
    bool fill = true;
    bool waveform = false;
    png.draw_curve(samples, crv.get_definition(), fill, waveform)
```

Lua : 
```Lua
    local crv = hc.hypercurve(2048, 0, {
      hc.segment(0.5, 1, hc.cubic())  
      hc.segment(0.5, 0, hc.cissoid(1))  
    })

    crv:ascii_display("MyHybridCurve", "half cubic, half cissoid", "*")
    local fill = true
    local is_waveform = false
    crv:write_as_png("path_to/curve.png", fill, is_waveform)
    crv:write_as_wav( "path/curve.wav" )

    -- rescale the curve
    crv:normalize_y(-1, 1)
    local samp = crv:get_sample_at(1024)
    local samps = crv:get_samples()
```

Csound : 
```Csound
    icrv = hc_hypercurve(2048, 0, 
        hc_segment(0.5, 1, hc_cubic_curve()), 
        hc_segment(0.5, 0, hc_cissoid_curve(1)))
    
    // Will run the curve in the time of the i event in an instrument
    kenv = hc_run(icrv, linseg:k(0, p3, 1))
```


## Import hypercurve

C++ : 
```c++
#include"hypercurve.h"
```
Lua : 
```Lua
-- In the following line, replace .so with your library extension (dll on Windows or dylib on Macos)
package.cpath = package.cpath .. ";/your/path/to/hypercurve/?.so"
local hc = require("liblua_hypercurve")
```
Csound : 
In csound you can manually import the library like below, or simply put the library in csound plugins path.
```Csound
<CsOptions>
--opcode-lib=/your/path/to/libcsound_hypercurve.so
</CsOptions>
```

## Hypercurve

C++ : 
```c++
    auto crv = hypercurve::hypercurve(size_t size_in_samples, double y_start, std::vector<hypercurve::segment> segment_list);    
    // With possible alias
    auto crv = hypercurve::curve(size_t size_in_samples, double y_start, std::vector<hypercurve::segment> segment_list);    
```
Lua : 
```Lua
    local crv = hc.hypercurve(integer size_in_samples, number y_start, table {segments})
```
Csound : 
```Csound
    icrv = hc_hypercurve(int isize_in_samples, float iy_start, isegment1 , [isegment2, isegment3...])
```

## Segment

C++ : 
```c++
    auto seg = hypercurve::segment(double fractional_size, double y_destination, std::shared_ptr<curve_base> curve_type);    
```
Lua : 
```Lua
    local seg = hc.segment(number fractional_size, number y_destination, curve_base)
```
Csound : 
```Csound
    iseg = hc_segment(float fractional_size, float y_destination, curve_base icrv_base)
```

<!---
C++ : 
```c++
    hypercurve::share( )
```
Lua : 
```Lua
    hc.
```
Csound : 
```Csound
    hc_
```
-->


## Curve Base

In Hypercurve, a Curve base represents the algorithm of a specific curve. Some of them take one or several constant parameters.

### Diocles cissoid curve

C++ : 
```c++
    hypercurve::share( hypercurve::diocles_curve(double a) );
    // Alias
    hypercurve::share( hypercurve::cissoid_curve(double a) );
```
Lua : 
```Lua
    hc.diocles(number a)
    -- Alias 
    hc.cissoid(a)
```
Csound : 
```Csound
    hc_diocles(float iarg_a)
    // Alias
    hc_cissoid(float iarg_a)
```

### Cubic curve

C++ : 
```c++
    hypercurve::share( hypercurve::cubic_curve() );
```
Lua : 
```Lua
    hc.cubic()
```
Csound : 
```Csound
    hc_cubic()
```

### Power curve

C++ : 
```c++
    hypercurve::share( hypercurve::power_curve(double power) );
```
Lua : 
```Lua
    hc.power(number power)
```
Csound : 
```Csound
    hc_power_curve(float ipower)
```
### Hamming / Hanning / Blackman curves

C++ : 
```c++
    hypercurve::share( hypercurve::hamming_curve() )
    hypercurve::share( hypercurve::hanning_curve() )
    hypercurve::share( hypercurve::blackman_curve() )
```
Lua : 
```Lua
    hc.hamming()
    hc.hanning()
    hc.blackman()
```
Csound : 
```Csound
    hc_hamming_curve()
    hc_hanning_curve()
    hc_blackman_curve()
```
### Gaussian curve (bell)

C++ : 
```c++
    hypercurve::share( hypercurve::gaussian_curve(double A, double c) );
    // Alias
    hypercurve::share( hypercurve::gauss_curve(double A, double c) );

```
Lua : 
```Lua
    hc.gaussian(number A, number c)
    -- Alias
    hc.gauss(number A, number c)
```
Csound : 
```Csound
    hc_gaussian_curve(float iA, float ic)
    // Alias
    hc_gauss_curve(float iA, float ic)
```

### Toxoid curve (duplicatrix cubic curve)

C++ : 
```c++
    hypercurve::share( hypercurve::toxoid_curve(double a) );
    // Alias
    hypercurve::share( hypercurve::duplicatrix_cubic_curve(double a) );
```
Lua : 
```Lua
    hc.toxoid(number a)
    -- Alias
    hc.duplicatrix_cubic(number a)
```
Csound : 
```Csound
    hc_toxoid_curve(float ia)
    // Alias
    hc_toxoid_curve(float ia)
```

### Catenary curve (funicular)

C++ : 
```c++
    hypercurve::share( hypercurve::catenary_curve(double a) );
    // Alias
    hypercurve::share( hypercurve::funicular_curve(double a) );
```
Lua : 
```Lua
    hc.catenary(number a)
    -- Alias
    hc.funicular(number a)
```
Csound : 
```Csound
    hc_catenary_curve(float ia)
    // Alias
    hc_funicular_curve(float ia)
```



### Tightrope Walker curve 

C++ : 
```c++
    hypercurve::share( hypercurve::tightrope_walker_curve(double a, double b) );
```
Lua : 
```Lua
    hc.tightrope_walker(number a, number b)
```
Csound : 
```Csound
    hc_tightrope_walker_curve(float ia, float ib)
```


### Quadratic Bezier curve
C++ : 
```c++
    hypercurve::share( hypercurve::quadratic_bezier_curve( hypercurve::control_point cp ) );
```
Lua : 
```Lua
    hc.quadratic_bezier( hc.control_point cp )
```
Csound : 
```Csound
    hc_quadratic_bezier_curve( hc_control_point cp )
```


### Cubic Bezier curve

C++ : 
```c++
    hypercurve::share( hypercurve::cubic_bezier_curve( hypercurve::control_point cp1, hypercurve::control_point cp2) );
```
Lua : 
```Lua
    hc.cubic_bezier(hc.control_point cp1, hc.control_point cp2)
```
Csound : 
```Csound
    hc_cubic_bezier_curve(hc_control_point cp1, hc_control_point cp2)
```


### Cubic spline curve

C++ : 
```c++
    hypercurve::share(hypercurve::cubic_spline_curve(std::vector<control_point> cp_list) );
```
Lua : 
```Lua
    hc.cubic_spline(table {hc.control_point})
```
Csound : 
```Csound
    // Not implemented yet
```


### Catmull Rom pline curve

C++ : 
```c++
    hypercurve::share( hypercurve::catmull_rom_spline_curve( hypercurve::control_point cp1,  hypercurve::control_point cp2) );
```
Lua : 
```Lua
    hc.catmull_rom_spline(hc.control_point cp1, hc.control_point cp2)
```
Csound : 
```Csound
    hc_catmull_rom_spline_curve(hc_control_point cp1, hc_control_point cp2)
```