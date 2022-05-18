
# Hypercurve for Faust  

You don't need to use cmake, since the whole library is contained as source only inside . The files `hypercurve.lib` and `hypercurve_faust.h` located 
inside [hypercurve_faust_lib](./hypercurve_faust_lib) folder can be copied to your Faust library directoy.

# Build 

If you still need to generate the output files yourself : 
The CMakeLists.txt file will generate a header only file that will later be included in the faustlibrary. In order to do that, you need to have [Quom](https://github.com/Viatorus/quom) installed.  
  
The Faust build basically consists of two steps :  
  
* Generate a header only hypercurve_faust.h that contains all hypercurve source code.  
* Move hypercurve.lib and hypercurve_faust.h to the build directory and to the output direcctory. 
  
Then, you should be able to copy these two files to your Faust library directory.

