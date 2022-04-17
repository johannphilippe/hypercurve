
# hypercurve


```
                                       awesome hybrid curve

	1 +--------------------------------------------------------------------------------+
          |        ********                                                                |
          |       *        **                                                              |
          |                                                                                |
          |      *           *******                                                       |
          |                         *                                                      |
          |     *                                     ******************                   |
          |                                      *****                  *****              |
          |                                    **                            ***           |
          |    *                     *      ***                                 ***        |
          |                               **                                       **      |
          |   *                          *                                           **    |
          |                             *                                              *   |
          |  *                         *                                                *  |
          | *                         *                                                  * |
          |*                                                                              *|
        0 +--------------------------------------------------------------------------------+
          0                                                                                 1

          +--------------------------------------------------------------------------------+
          |   * A combination of gaussian, cissoid, power of 9, and cubic bezier curve     |
          +--------------------------------------------------------------------------------+

```




# What is it ? 


A set of classes to create hybrid curves and use them for audio and music (or whatever). 
It can be used in C++, Lua or Csound.

On a mathematical point of view, it has to be noticed that since curves needs to be scaled in and mapped together, there must be a few approximations inherent to the application.

## Implemented curves 


- Cissoid (Diocles curve) 
- Cubic 
- Power curve (choose your power of x)
- Bezier (Cubic & Quadratic)
- Cubic Spline - Not implemented in Csound yet
- Catmull Rom Spline
- Hanning / Hamming / Blackamn
- Gauss 
- Toxoid (aka duplicatrix_cubic)
- Catenary (aka funicular)
- Tightrope Walker curve  

- Typed curves : inspired from Csounds [GEN16](http://www.csounds.com/manual/html/GEN16.html)
- User defined curves - pass it a function (or a lambda in C++), that returns y for any x between 0 and 1. Not implemented in Csound.


## How to use it 


There are three ways to use it : in C++, Csound or in Lua. Cmake will help you build libraries that can be used in those languages. You will find C++ example under hypercurve_test/test.cpp, Csound example under csound_opcode/test.csd, and Lua example under lua_module/test.lua. 

The basic syntax stands as follow : 
* `hypercurve(integer size, double y_start, {segment_list});`
Where `size` is the size in samples, `y_start` is the starting point of the curve, and segment list is a list of segments. 
*  `segment(double frac, double y_destination, curve_type crv);`
Where `frac` is the fractional size of the segment (fraction between 0 and 1), `y_destination` is the target point, and `crv`  a curve picked from hypercurve.


## A simple C++ example 

```c++
#include"hypercurve.h"

using namespace hypercurve;
const int definition = 16384;
double y_start = 0;
curve c(definition, y_start, 
	{
		// segment(fractional_size, y_destination, curve
		segment(frac(1,2), 1.0, share(cissoid_curve(1))),
		segment(0.5, 0.0, share(blackman_curve()))
	}); 

// Then access samples with double *get_samples() 
c.get_samples();
```

## A simple Csound example

```csound
instr 1
	icrv = hypercurve(2048, 0, 
				segment(1/2, 1, diocles_curve(1)),
				segment(1/2, 0, hanning_curve()))
	kenv = run_hypercurve(icrv, linseg(0, p3, 1))
	ao = vco2(0.3, 300) * kenv
	outs(ao, ao)
endin

```

## A simple Lua example

```lua
package.cpath = package.cpath .. ";/your/path/to/hypercurve/?.so;"
local hc = require("liblua_hypercurve")

local definition = 16384
local y_start = 0
local crv = hc.hypercurve(definition, y_start, 
	{
		hc.segment(1/2, 1.0, hc.cissoid(1.0)),
		hc.sement(1/2, 0.0, hc.cubic(0.0))
	})

// Write as 24 bits 48KHz wav
hc.write_as_wav("path/to/outfile.wav", crv)
```

## Features 

The  `curve`  takes a list of segments, each having a fractional size. If the sum of all these segments is not exactly one, they will be rescaled so that they fit 0-1 range. 




# Build


First clone the repo with submodules : 
``` git clone https://github.com/johannphilippe/hypercurve.git --recurse-submodules ```

You should check that Lua is installed on your system. If it is not, or if compilation retrns error, you should install a Lua 5.1 version to the standard installation path. Make sure you have the dynamic library installed, and the headers `lauxlib.h` and `lua.h` are available on your system.

Then : 
```
cd hypercurve
mkdir build && cd build
cmake .. -DBUILD_ALL=TRUE
make
```

If you just want to build Lua module or Csound opcode, then just use 

```
cmake .. -DBUILD_CSOUND_OPCODE=TRUE
cmake .. -DBUILD_LUA_MODULE=TRUE
```

On some platforms (e.g. Windows) you might need to set the Lua paths with the following options :
```
cmake .. -DBUILD_LUA_MODULE=TRUE -DLUA_INCLUDE_DIR=/you/dir/include -DLUA_LIBRARIES=/path/to/lua.lib
```


The PNG writer [fpng](https://github.com/richgel999/fpng) used for hypercurve has SSE support. This can be enabled with `-DSSE=1`.



# TODO

* REAPER/Reascript -> see https://forum.cockos.com/showthread.php?p=2543755#post2543755

* Implement a polynomial with varargs (like a, b, c  : ax^3 + bx^2 + c)

* A real good picture in README to show what it actually looks like

* Operator for curves in Lua and Csound. In Csound it should be functions (sum_hypercurve, substract, multiply, div) that create a new one. In Lua just a layer on top of C++ operators with metatables.
* Csound RT opcode : cubic spline (mem alloc). Also missing X rescale, y_rescale, and user_defined curve
* Implement rescale (Lua, Csound) a

* Lagrange interpolation for curve extraction.

* Documentation : add pictures of each curve.

* Hard one -> all curves allowing one sample processing (including cubic spline) to allow no-table processing.

* Semantic : keep curve in "diocles_curve" or remove it ? (choose btw lua or csound style)

## Curves to implement

* Cardioid / hypercardioid
* Elastic curve : https://mathcurve.com/courbes2d.gb/linteaire/linteaire.shtml
* Simple log/exp ?
* Kulp quartic

* Ideas here https://mathcurve.com/courbes2d/courbes2d.shtml

# External libraries

This library uses libsndfile as an external submodule.
Currently, [libsndfile](https://github.com/libsndfile/libsndfile) is only used in Lua and in test.cpp. So hypercurve C++ library itself does not have any external dependency.

It also includes source files from several open-source projects : 
*  [AsciiPlot](https://github.com/joehood/asciiplotter) source code with license under src/asciiplot folder.
* [lua-compat-5.3](https://github.com/keplerproject/lua-compat-5.3) which provides an API compatibility from 5.1 to 5.3
* [fpng](https://github.com/richgel999/fpng) - a great C++ PNG reader/writer.

