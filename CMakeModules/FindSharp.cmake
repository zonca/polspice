#
#
#  CMake module looking for sharp library 
#  (compiled separately in HEALPix 3.60)
#  required when linking PolSpice with HEALPix
#  2019-06-11: v1.0, E. Hivon
#
find_library(SHARP_LIBRARY NAMES sharp)
set(SHARP_LIBRARIES ${SHARP_LIBRARY})
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SHARP DEFAULT_MSG SHARP_LIBRARIES)

mark_as_advanced(
SHARP_LIBRARY
SHARP_LIBRARIES
)
