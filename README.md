# hypercurve

A set of polynomial curves functions that can be used as audio control curves. 

# Build

First clone the repo with submodules : 
``` git clone repo --recurse-submodules ```

You should check that Lua is installed on your system. 
Then : 
```
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
* Bezier cubic / quadratic
* Spline (tricky, but feasable)


# External libraries

This library uses libsndfile as an external submodule.
