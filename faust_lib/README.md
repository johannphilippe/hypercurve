
# Hypercurve for Faust  
  
The CMakeLists.txt file will generate a header only file that will later be included in the faustlibrary. In order to do that, you need to have [Quom](https://github.com/Viatorus/quom) installed.  
  
The Faust build basically consists of two steps :  
  
* Generate a header only hypercurve_faust.h that contains all hypercurve.  
* Move hypercurve.lib and hypercurve_faust.h to the build directory.  
  
Then, you should be able to copy these two files to your Faust library directory.

