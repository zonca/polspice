#
# CMake file for PolSpice
# To use it do
# mkdir build ; cd build ; cmake ..  ; make
#
# tested successfully with versions 3.8.1, 3.9.6, 3.10.* and 3.15.3 of cmake
# cmake 3.7.1 only works with libcfitsio.a (ie, static, NOT shared/dynamic)
#
cmake_minimum_required(VERSION 3.7)
project(spice)
set(VERSION "3.7.0")

include(CheckIncludeFile)
include(CheckLibraryExists)
include(FindPkgConfig)
include(FindPackageHandleStandardArgs)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)

if (NOT HEALPIX)
    set (HEALPIX "$ENV{HEALPIX}")
endif()

get_filename_component (C_COMPILER_NAME ${CMAKE_C_COMPILER} NAME)
enable_language (Fortran)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

message("  ================================
  CMakeLists.txt for Spice ${VERSION}.")
message("  Assuming HEALPix library (>=3.60) to be in ${HEALPIX}/lib${SUFFIX}/*
  F90/95 compiler: ${CMAKE_Fortran_COMPILER}
  C      compiler: ${CMAKE_C_COMPILER}
  For more information see INSTALL")
message("  ================================")

# set(HEALPIX "" CACHE PATH "spice: Specify the HEALPix directory location (overrules the env. variable $HEALPIX)")
# set(HEALPIX_LIB "" CACHE PATH "spice: location of HEALPix library (overrules $HEALPIX/lib$SUFFIX)")
# set(HEALPIX_INCLUDE "" CACHE PATH "spice: location of HEALPix modules (overrules $HEALPIX/include$SUFFIX)")
# set(SUFFIX "" CACHE PATH "spice: libraries will be looked for in $HEALPIX/lib$SUFFIX and modules in $HEALPIX/include$SUFFIX")
# set(SHARPDIR "" CACHE PATH "spice: location of libsharp library needed for HEALPix >= 3.60 (overrules $HEALPIX/lib)")
# set(CFITSIO "" CACHE PATH "spice: location of cfitsio library (can be a directory or the full path)")

# FFLAGS depend on the compiler
if (Fortran_COMPILER_NAME MATCHES "gfortran.*")  # gfortran (GNU)
  set (CMAKE_Fortran_FLAGS "-fopenmp -O3")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*") # ifort (intel)
  set (CMAKE_Fortran_FLAGS "-qopenmp -O3")
  #set (CMAKE_Fortran_FLAGS "-qopenmp -O3 -check all -traceback -ftrapuv")
elseif (Fortran_COMPILER_NAME MATCHES "g95.*")   # g95 (no OpenMP support)
  set (CMAKE_Fortran_FLAGS " -w -ffree-form -O3")
elseif (Fortran_COMPILER_NAME MATCHES "pgf9.*")  # pgf90, pgf95 (PGI, formerly Portland Group)
  set (CMAKE_Fortran_FLAGS " -mp -O3")
elseif (Fortran_COMPILER_NAME MATCHES "crayf.*") # crayftn (Cray)
  set (CMAKE_Fortran_FLAGS " -h omp -O3")
else ()
  message ("
    CMAKE_Fortran_COMPILER full path: ${CMAKE_Fortran_COMPILER}
    Fortran compiler: ${Fortran_COMPILER_NAME}
     Use the environment variable   FC
    or the cmake flag   -DCMAKE_Fortran_COMPILER=xxx
    to specify the compiler")
  message (FATAL_ERROR "Compiler not known")
endif ()

# Option -DDOUBLE=1
option (DOUBLE "use double precision variables for maps and alm" ON)
if (DOUBLE)
  string(APPEND CMAKE_Fortran_FLAGS " -DDOUBLE=1 ")
endif ()

message(" * Compiling with:   (${Fortran_COMPILER_NAME}) ${CMAKE_Fortran_COMPILER} ${CMAKE_Fortran_FLAGS}")
message(" *            and:   (${C_COMPILER_NAME}) ${CMAKE_C_COMPILER} ${CMAKE_C_FLAGS}")


# define directories and executables
set(BIN_DIR ${PROJECT_SOURCE_DIR}/bin)
set(SRC_DIR ${PROJECT_SOURCE_DIR}/src)
if (NOT EXEC)
   set(EXEC spice)
endif()
set(BLD_DIR ${PROJECT_BINARY_DIR})
message (" * source directory: ${SRC_DIR} ")
message (" * executable:       ${EXEC} ")

# define target
set (TARGET ${EXEC})


# define source files for target
set(SRCS
  ${SRC_DIR}/spice_parameters.F90
  ${SRC_DIR}/spice_common.F90
  ${SRC_DIR}/fits_spice.F90
  ${SRC_DIR}/apodize_mod.f90 ${SRC_DIR}/convert_to_ascii.f90 ${SRC_DIR}/do_legendre.f90 
  ${SRC_DIR}/rec3jj.f90      ${SRC_DIR}/rec3jjcmb.F90        ${SRC_DIR}/remove_dipole_mod.f90
  ${SRC_DIR}/compute_windows.F90
  ${SRC_DIR}/deal_with_xi_and_cl.f90
  ${SRC_DIR}/cumul2.f90
  ${SRC_DIR}/spice_subs.f90
  ${SRC_DIR}/deal_with_filenames.F90 ${SRC_DIR}/deal_with_files.F90 ${SRC_DIR}/deal_with_options.F90
  ${SRC_DIR}/spice.F90
    ${SRC_DIR}/sets.F90 ${SRC_DIR}/map_type.F90 ${SRC_DIR}/alm_type.F90 
    ${SRC_DIR}/sharp_f_interface.F90 ${SRC_DIR}/sharp_interface.c )

# create target  
add_executable(${TARGET} ${SRCS})

# deal with SUFFIX
# if (SUFFIX)
#     set (SUFF "${SUFFIX}")
# else ()
#     set (SUFF "")
# endif()
###message(" >>>>>>        ${SUFF}   ${SUFFIX}  ")


# link with HEALPix lib
find_library(HEALPIX_LIB healpix PATHS ${HEALPIX_LIB}   NO_DEFAULT_PATH)
find_library(HEALPIX_LIB healpix PATHS "${HEALPIX}/lib${SUFFIX}" NO_DEFAULT_PATH)
#find_library(HEALPIX_LIB healpix PATHS "$ENV{HEALPIX}/lib${SUFF}" NO_DEFAULT_PATH)
find_library(HEALPIX_LIB healpix)
if (HEALPIX_LIB)
  message(" * Linking to ${HEALPIX_LIB}")
  add_library(healpix UNKNOWN IMPORTED)
  target_link_libraries (${TARGET} ${HEALPIX_LIB})
else()
  message(FATAL_ERROR "
  HEALPix library not found in ${HEALPIX}/lib${SUFFIX}/*
  Use the  cmake flags  
   -DHEALPIX=/correct/path/to/healpix/top/directory  (and -DSUFFIX=suffix)
    to use {HEALPIX}/lib{SUFFIX}/*   instead
   or 
   -DHEALPIX_LIB=/correct/path/to/healpix/library/directory 
    to use {HEALPIX_LIB}/*
  For more information see INSTALL")
endif()

# 'include' HEALPix modules
if (NOT HEALPIX_INCLUDE)
  set(HEALPIX_INCLUDE "${HEALPIX}/include${SUFFIX}" )
#   if (HEALPIX)
#     set(HEALPIX_INCLUDE "${HEALPIX}/include${SUFF}" )
#   else()
#     set(HEALPIX_INCLUDE "$ENV{HEALPIX}/include${SUFF}" )
#   endif()
  set(HEALPIX_C_INCLUDE "${HEALPIX}/include" )
endif()
message(" * Include files from ${HEALPIX_INCLUDE}")
target_include_directories (${TARGET} PRIVATE ${HEALPIX_INCLUDE})
target_include_directories (${TARGET} PRIVATE ${HEALPIX_C_INCLUDE})

# link with libsharp:
#   look in ${HEALPIX}/lib or in user provided SHARPDIR
if (NOT SHARPDIR)
  set(SHARPDIR "${HEALPIX}/lib" )
endif()
set(SHARPINC "${SHARPDIR}/../include" )
set(CMAKE_PREFIX_PATH ${SHARPDIR})
find_package(Sharp QUIET)# use ${CMAKE_MODULE_PATH}/FindSharp.cmake
if (SHARP_FOUND)
    target_link_libraries (${TARGET} ${SHARP_LIBRARIES})
    message(" * Linking to ${SHARP_LIBRARIES}")
    find_path(SHARP_INCLUDES
	  NAMES c_utils.h sharp.h sharp_geomhelpers.h sharp_almhelpers.h
	  PATHS ${SHARPINC} ${SHARPINC}/libsharp ${SRC_DIR}
	  PATH_SUFFIXES NO_DEFAULT_PATH)
    target_include_directories (${TARGET} PRIVATE ${SHARP_INCLUDES})
    message(" * Include files from ${SHARP_INCLUDES}")
endif()
#	  NAMES sharp.h sharp_geomhelpers.h sharp_almhelpers.h c_utils.h

# link with cfitsio
if (CFITSIO)
  #message ("User Cfitsio ${CFITSIO}")
  if (${CFITSIO} MATCHES "(.*\.a|.*\.so|.*\.dylib)$")
    set (CFITSIO_LIBRARY "${CFITSIO}") # user specified file
  else()
    find_library(CFITSIO_LIBRARY cfitsio PATHS ${CFITSIO} NO_DEFAULT_PATH) # user specified directory
  endif()
else()
  #message ("Default Cfitsio")
  find_library(CFITSIO_LIBRARY cfitsio) # default locations
endif()
if (CFITSIO_LIBRARY)
  message(" * Linking to ${CFITSIO_LIBRARY}")
  target_link_libraries (${TARGET} ${CFITSIO_LIBRARY})
else()
  message(FATAL_ERROR "
  The cfitsio library was not found in usual places
  Use -DCFITSIO=cfitsio_directory to specify the location of the cfitsio library.")
endif()
# recent versions of cfitsio use curl if it is available


if (CURL)
  if (${CURL} MATCHES "(.*\.a|.*\.so|.*\.dylib)$")
    set (CURL_LIBRARY "${CURL}") # user specified file
  else()
    find_library(CURL_LIBRARY curl PATHS ${CURL} NO_DEFAULT_PATH) # user specified directory
  endif()
  if (CURL_LIBRARY)
     target_link_libraries (${TARGET} ${CURL_LIBRARY})
     message(" * Linking to ${CURL_LIBRARY}")
  else()
     message("curl not found.")
  endif()
else()
  find_package(CURL)
  if (CURL_FOUND)
     target_link_libraries (${TARGET} curl)
     message(" * Linking to ${CURL_LIBRARIES}")
  else()
     message("curl not found.
  Use -DCURL=curl_directory to specify the location of the curl library.")
  endif()
endif()


file(GLOB MY_PRO  "${SRC_DIR}/*.pro")
file(GLOB MY_PY   "${SRC_DIR}/*.py")
file(COPY ${MY_PRO} ${MY_PY} ${MY_CODE} DESTINATION ${BIN_DIR})

add_custom_command(TARGET ${TARGET}
          POST_BUILD
          COMMAND ${CMAKE_COMMAND} -E make_directory ${BIN_DIR}
          COMMAND ${CMAKE_COMMAND} -E copy ${TARGET} ${BIN_DIR}
)


#----------- targets for cleaning ----------

# Create a 'make vclean'
# removing the intermediate files and the cmake specific files from build directory
# Does not remove the whole build directory in case build=src
add_custom_target(vclean
  COMMAND ${CMAKE_MAKE_PROGRAM} clean
  COMMAND ${CMAKE_COMMAND} -E remove      ${BLD_DIR}/cmake_install.cmake 
					  ${BLD_DIR}/Makefile 
					  ${BLD_DIR}/CMakeCache.txt 
  COMMAND ${CMAKE_COMMAND} -E remove_directory      ${BLD_DIR}/CMakeFiles
)

# Create a 'make tidy'
# doing  a 'make vclean' and removing the bin directory
add_custom_target(tidy
  COMMAND ${CMAKE_COMMAND} -E remove_directory      ${BIN_DIR}
  COMMAND ${CMAKE_MAKE_PROGRAM} vclean
)




