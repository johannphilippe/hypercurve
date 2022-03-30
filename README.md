# hypercurve

A set of polynomial curves functions that can be used as audio control curves. 



# What is it ? 


A set of classes to create composite curves and export them in various format. 
Currently, you can export curves as `.wav` files. 

## Implemented curves 


- Cissoid (Diocles) 
- Cubic 
- Bezier (Cubic & Quadratic)



## How to use it 


There are two ways to use it : in C++ or in Lua. The CMake build system builds two libraries that can be used in both languages. You will find C++ example under hypercurve_test/test.cpp, and Lua example under lua_module/test.lua.  



# Build


First clone the repo with submodules : 
``` git clone https://github.com/johannphilippe/hypercurve.git --recurse-submodules ```

You should check that Lua is installed on your system. If it is not, or if compilation returns error, you should install a Lua 5.1 version to the standard installation path. Make sure you have the dynamic library installed, and the headers `lauxlib.h` and `lua.h` are available on your system.

Then : 
```
cd hypercurve
mkdir build && cd build
cmake ..
make
```


# TODO

* Ruby gem config
* Lua Improve (OOP for curve class) and tests (Reaper)
* Csound opcode config

## TODO polynomials

* Cardioid / hypercardioid
* Spline (tricky, but feasable)


# External libraries

This library uses libsndfile as an external submodule.
