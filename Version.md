# 0.0.1 

First release

# 0.0.2 

## New 

- Faust (FFI) implementation
- "_curve" suffix for all curve algorithm (curve_base) in every frontend (including Lua)
- Cubic spline curve ported to Csound opcodes
- Added logarithmic_curve and exponential_curve approximations
- Renamed `normalize_y` to `scale`. Added `normalize` and `norm` (aliases for `scale(0,1)`).
- Added mathematic operators to operate curve with numbers, where second argument is the number 
  - in Csound and Faust it is `hc_addn(curve, 2)` to add a number to a curve. 
  - In Lua and C++, simply use operators `curve + 2`. 
- Removed Sndfile dependency

## Fixed

- Cubic spline implementation improved
- Renaming mouse curve to "mouth curve" - (typo fix)
- Round instead of floor for Csound index calculations
- Fixed size issues with 0 padding when one sample is missed (fractional size divided by 3, 6, 9 ...)
- Fixed Cubic Spline allocation issues
