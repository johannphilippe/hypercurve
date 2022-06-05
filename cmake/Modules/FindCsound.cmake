# Try to find the Csound library.
# Once done this will define:
#  CSOUND_FOUND - System has the Csound library
#  CSOUND_INCLUDE_DIRS - The Csound include directories.
#  CSOUND_LIBRARIES - The libraries needed to use the Csound library.

if(APPLE)
    find_path(CSOUND_INCLUDE_DIR csound.h HINTS /Library/Frameworks/CsoundLib64.framework/Headers
        "$ENV{HOME}/Library/Frameworks/CsoundLib64.framework/Headers")
else()
    find_path(CSOUND_INCLUDE_DIR csound.h PATH_SUFFIXES csound)
endif()

if(APPLE)
    find_library(CSOUND_LIBRARY NAMES CsoundLib64 HINTS /Library/Frameworks/CsoundLib64.framework/
        "$ENV{HOME}/Library/Frameworks/CsoundLib64.framework")
else()

    find_library(CSOUND_LIBRARY NAMES csound64 csound)
endif()

if(APPLE)
    find_path(CSOUND_PLUGINS_DIR libarrayops.dylib HINTS /Library/Frameworks/CsoundLib64.framework
        "$ENV{HOME}/Library/Frameworks/CsoundLib64.framework")
elseif(WIN32)
    find_path(CSOUND_PLUGINS_DIR arrayops.dll HINTS "C:/Program Files/Csound6_x64" "C:/Program Files/Csound6_x64/plugins64")
else()
    find_path(CSOUND_PLUGINS_DIR libarrayops.so HINTS /usr/local/lib/csound /usr/lib/csound)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set CSOUND_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CSOUND
                                  CSOUND_LIBRARY CSOUND_INCLUDE_DIR CSOUND_PLUGINS_DIR)
mark_as_advanced(CSOUND_INCLUDE_DIR CSOUND_LIBRARY)

set(CSOUND_INCLUDE_DIRS ${CSOUND_INCLUDE_DIR})
set(CSOUND_LIBRARIES ${CSOUND_LIBRARY} )
