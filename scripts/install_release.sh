#!/bin/bash

# Not ready yet

case "$OSTYPE" in
  solaris*)
	  echo "SOLARIS" ;;
  darwin*) 
	  archive="hypercurve_MACOS_Clang12.zip"
	  wget "https://github.com/johannphilippe/hypercurve/releases/latest/download/hypercurve_MACOS_Clang12.zip"
	  echo "OSX" ;; 
  linux*)
	  archive="hypercurve_Linux_Clang10.zip"
	  wget "https://github.com/johannphilippe/hypercurve/releases/latest/download/hypercurve_Linux_Clang10.zip"
       	  echo "LINUX" ;;
  bsd*)    
	  echo "BSD" ;;
  msys*)   
	  archive="hypercurve_Windows11_MSVC19_x64.zip"
	  wget "https://github.com/johannphilippe/hypercurve/releases/latest/download/hypercurve_Windows11_MSVC19_x64.zip"
	  echo "WINDOWS" ;;
  cygwin*) 
	  echo "ALSO WINDOWS" ;;
  *)       
	  echo "unknown: $OSTYPE" ;;
esac

unzip $archive
