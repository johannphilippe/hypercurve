
# hypercurve


```
                                            CatmullRom

        1 +--------------------------------------------------------------------------------+
          |                                  ********                                      |
          |                                **        ***                                   |
          |                             ***             ***                                |
          |                           **                   *****                           |
          |                         **                          *******                    |
          |                        *                                   ****                |
          |                      **                                        ***             |
          |                    **                                             **           |
          |                  **                                                 **         |
          |                **                                                     **       |
          |              **                                                         **     |
          |           ***                                                             *    |
          |        ***                                                                 **  |
          |     ***                                                                      * |
          |*****                                                                          *|
        0 +--------------------------------------------------------------------------------+
          0                                                                                 1

          +--------------------------------------------------------------------------------+
          |   * y = catmullrom(X)                                                          |
          +--------------------------------------------------------------------------------+

```


# What is it ? 


A set of classes to create hypbrid curves and use them for audio and music (or whatever). 
Currently, you can export curves as `.wav` files. 

## Implemented curves 


- Cissoid (Diocles) 
- Cubic 
- Bezier (Cubic & Quadratic)
- Cubic Spline
- Catmull Rom Spline

## How to use it 


There are three ways to use it : in C++, Csound or in Lua. Cmake will help you build libraries that can be used in those languages. You will find C++ example under hypercurve_test/test.cpp, Csound example under csound_opcode/test.csd, and Lua example under lua_module/test.lua. 


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
		segment(0.5, 0.0, share(cubic_curve()))
	}); 

// Then access samples with double *get_samples() 
c.get_samples();
```

## A simple Csound example

```csound
instr 1
	icrv = hypercurve(2048, 0, 
				segment(1/2, 1, diocles_curve(1)),
				segment(1/2, 0, diocles_curve(1)))
	kenv = run_hypercurve(icrv, linseg(0, p3, 1))
	ao = vco2(0.3, 300) * kenv
	outs(ao, ao)
endin

```

## A simple Lua example

```lua
package.cpath = package.cpath .. ";/your/path/to/hypercurve"
local hc = require("liblua_hypercurve")

local definition = 16384
local y_start = 0
local crv = hc.curve(definition, y_start, 
	{
		hc.segment(1/2, 1.0, hc.cissoid_curve(1.0)),
		hc.segment(1/2, 0.0, hc.cubic_curve(0.0))
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
cmake ..
make
```

# TODO

* REAPER/Reascript -> see https://forum.cockos.com/showthread.php?p=2543755#post2543755

* Csound RT opcode : cubic spline (mem alloc). Also missing X rescale, y_rescale.

* Lagrange interpolation for curve extraction.

* Once all done -> a proper documentation for C++/Lua, and Csound implementations

## Curves to implement

* Cardioid / hypercardioid
* Hanning, Hamming, Blackmann (half)


# External libraries

This library uses libsndfile as an external submodule.
It also includes source files from several open-source projects : 
*  [AsciiPlot](https://github.com/joehood/asciiplotter) source code with license under src/asciiplot folder.
* [lua-compat-5.3](https://github.com/keplerproject/lua-compat-5.3) which provides an API compatibility from 5.1 to 5.3

