
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

Hypercurve is a library allowing you to use several curve algorithms into a single 2D envelope. It is designed to be used in audio applications, for people who know how to enjoy a finely shaped curve. 
As shown above, you can perfectly combine gaussian curve with bezier curve, and plenty of other curve algorithms. 
It can be used in C++, Lua or Csound.

Every curve algorithm is different. In an audio application, you can assign an envelope to any kind of parameter. For some parameters (like frequency), the way a value goes up and down in time changes everything to the perception of sound. Thus, the possibility of having a finely shaped envelope is truly essential. This, is the purpose of Hypercurve. 

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
	icrv = hc_hypercurve(2048, 0, 
				hc_segment(1/2, 1, hc_diocles_curve(1)),
				hc_segment(1/2, 0, hc_hanning_curve()))
	kenv = hc_run(icrv, linseg(0, p3, 1))
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
* Fix toxoid curve (hybrid example)
* Fix bezier and spline control points (is it relative or absolute to the segment ? think about catmull rom)
* Lua semantics : append "curve"  to curve_base methods
* REAPER/Reascript -> see https://forum.cockos.com/showthread.php?p=2543755#post2543755
* A real good picture in README to show what it actually looks like
* Lagrange interpolation for curve extraction.
* Documentation : add pictures of each curve.
* Hard one -> all curves allowing one sample processing (including cubic spline) to allow no-table processing.
## Curves to implement
* Cardioid / hypercardioid
* Elastic curve : https://mathcurve.com/courbes2d.gb/linteaire/linteaire.shtml
* Simple log/exp ?
* Kulp quartic
* Puntiforme https://mathcurve.com/courbes2d/puntiforme/puntiforme.shtml
* Mouse https://mathcurve.com/courbes2d/bouche/bouche.shtml
* Bicorn AKA cocked hat  https://mathcurve.com/courbes2d/bicorne/bicorne.shtml
* Legendre polynome
* Ideas here https://mathcurve.com/courbes2d/courbes2d.shtml
# External libraries
This library uses libsndfile as an external submodule.
Currently, [libsndfile](https://github.com/libsndfile/libsndfile) is only used in Lua and in test.cpp. So hypercurve C++ library itself does not have any external dependency.
It also includes source files from several open-source projects : 
*  [AsciiPlot](https://github.com/joehood/asciiplotter) source code with license under src/asciiplot folder.
* [lua-compat-5.3](https://github.com/keplerproject/lua-compat-5.3) which provides an API compatibility from 5.1 to 5.3
* [fpng](https://github.com/richgel999/fpng) - a great C++ PNG reader/writer.``
