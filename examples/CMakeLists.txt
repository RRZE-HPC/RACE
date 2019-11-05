cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

# project name
project(RACE_example CXX C)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_PREFIX_PATH ${RACE_DIR})
# find RACE
find_package(RACE REQUIRED CONFIG)
message(STATUS "Found RACE")

#likwid library
link_directories($ENV{LIKWID_LIBDIR})
link_libraries(liblikwid.so)
include_directories($ENV{LIKWID_INCDIR})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DLIKWID_PERFMON -DLIKWID_MEASURE")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLIKWID_PERFMON -DLIKWID_MEASURE")


set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${RACE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${RACE_CXX_FLAGS}")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -xHOST -save-temps -fsource-asm")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xHOST -save-temps -fsource-asm")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -xCORE-AVX2  -save-temps -fsource-asm")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xCORE-AVX2  -save-temps -fsource-asm")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wno-write-strings")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-write-strings")
endif()

find_package(OpenMP)

# TODO: would be nicer if the race config.cmake file could provide te OpenMP flags
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(RACE_USE_SPMP 0 CACHE BOOL "Use RCM by SpMP for matrix preprocessing")

if (RACE_USE_SPMP)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(RACE_HAVE_SPMP 1)
        set(SPMP_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
        set_property(DIRECTORY PROPERTY EP_STEP_TARGETS build)
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            ExternalProject_Add(SpMP
                PREFIX ${PROJECT_BINARY_DIR}/SpMP
                URL https://github.com/IntelLabs/SpMP/archive/master.zip
                CONFIGURE_COMMAND ""
                BUILD_COMMAND make 
                BUILD_IN_SOURCE 1
                INSTALL_COMMAND ""
                )
        else() # GNU or Clang: Set CUSTOM_CXX variable in SpMP make
            ExternalProject_Add(SpMP
                PREFIX ${PROJECT_BINARY_DIR}/SpMP
                URL https://github.com/IntelLabs/SpMP/archive/master.zip
                CONFIGURE_COMMAND ""
                BUILD_COMMAND make CUSTOM_CXX=${CMAKE_CXX_COMPILER} 
                BUILD_IN_SOURCE 1
                INSTALL_COMMAND ""
                )
        endif()
    else()
        message(FATAL_ERROR "SpMP will not build with a compiler which is not Clang/GNU/Intel. Please disable SpMP or use a different compiler!")
    endif()
endif()


set(CBLAS_INCLUDE_DIR CBLAS_INCLUDE_DIR-NOTFOUND CACHE PATH "Location of *cblas.h")
# firstly, try to find a Cblas header in the default locations
if (NOT CBLAS_INCLUDE_DIR)
    message(STATUS "Try to find MKL in default paths and via MKLROOT")
    find_path(CBLAS_INCLUDE_DIR mkl_cblas.h $ENV{MKLROOT}/include)
endif()
if (NOT CBLAS_INCLUDE_DIR)
    #Have to search GSL with find_package because it will not be found by find_path
    message(STATUS "Try to find GSL in default paths and via find_package")
    find_package(GSL)
    find_path(CBLAS_INCLUDE_DIR gsl_cblas.h ${GSL_INCLUDE_DIRS}/gsl)
endif()
if (NOT CBLAS_INCLUDE_DIR)
    message(STATUS "Try to find default Cblas")
    find_path(CBLAS_INCLUDE_DIR cblas.h)
endif()

if (NOT CBLAS_INCLUDE_DIR)
    message(FATAL_ERROR "CBLAS_INCLUDE_DIR not set!")
endif()

include_directories(SYSTEM ${CBLAS_INCLUDE_DIR})

find_file(CBLAS_H_FULL NAMES "mkl_cblas.h" "cblas.h"  "gsl_cblas.h" CBLAS_H_FULL PATHS ${CBLAS_INCLUDE_DIR} PATH_SUFFIXES "cblas")
get_filename_component(CBLAS_H ${CBLAS_H_FULL} NAME)

if (${CBLAS_H} STREQUAL "mkl_cblas.h")
    message(STATUS "MKL BLAS found")
    set(RACE_HAVE_MKL 1)
    get_filename_component(MKL_LIBDIR ${CBLAS_INCLUDE_DIR}/../lib/ REALPATH)
    if (RACE_MIC_CROSSCOMPILE)
        set(MKL_PATH_SUFFIXES "mic")
    else ()
        set(MKL_PATH_SUFFIXES "intel64")
    endif()

    set(LIBMKL_CORE LIBMKL_CORE-NOTFOUND)
    set(LIBMKL_X_THREAD LIBMKL_X_THREAD-NOTFOUND)
    set(LIBMKL_INTEL_XLP64 LIBMKL_INTEL_XLP64-NOTFOUND)

    find_library(LIBMKL_CORE mkl_core PATHS ${MKL_LIBDIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})

    if (RACE_IDX64_LOCAL)
        add_definitions(-DMKL_ILP64)
        find_library(LIBMKL_INTEL_XLP64 mkl_intel_ilp64 PATHS ${MKL_LIBDIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
    else()
        find_library(LIBMKL_INTEL_XLP64 mkl_intel_lp64 PATHS ${MKL_LIBDIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
    endif()

    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel") 
        find_library(LIBMKL_X_THREAD mkl_intel_thread PATHS ${MKL_LIBDIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") 
        find_library(LIBMKL_X_THREAD mkl_gnu_thread PATHS ${MKL_LIBDIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
    else ()
        find_library(LIBMKL_X_THREAD mkl_gnu_thread PATHS ${MKL_LIBDIR} PATH_SUFFIXES ${MKL_PATH_SUFFIXES})
        #            message(FATAL_ERROR "OpenMP only possible with Intel or Gnu!")
    endif()

    set (BLAS_LIBRARIES "${LIBMKL_INTEL_XLP64};${LIBMKL_CORE};${LIBMKL_X_THREAD}")
    mark_as_advanced(LIBMKL_CORE LIBMKL_INTEL_XLP64 LIBMKL_X_THREAD) 
elseif (${CBLAS_H} STREQUAL "gsl_cblas.h")
    set(RACE_HAVE_GSL 1)
    get_filename_component(GSL-CONFIG_PATH ${CBLAS_INCLUDE_DIR}/../../bin REALPATH)
    list(APPEND CMAKE_PROGRAM_PATH ${GSL-CONFIG_PATH})
    find_package(GSL REQUIRED)
    message(STATUS "GSL BLAS found")
    set (BLAS_LIBRARIES "${GSL_LIBRARIES}" CACHE STRING "BLAS libraries")
else () #generic BLAS
    set (BLAS_LIBRARIES "" CACHE STRING "BLAS libraries")
    message(STATUS "Generic BLAS. BLAS_LIBRARIES has to be set manually before linking.")
endif()


configure_file (
    "${PROJECT_SOURCE_DIR}/include/config_eg.h.in"
    "${PROJECT_BINARY_DIR}/include/config_eg.h"
    )

#include_directories(${RACE_INCLUDE_DIR})
link_libraries(${RACE_LIBRARIES})
set(CXX_STANDARD 11)
Set(CXX_STANDARD_REQUIRED On)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${RACE_INCLUDE_DIRS})
include_directories(${PROJECT_BINARY_DIR}/include)

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/helpers/*.cpp)
add_executable(race ${SOURCES}  ${PROJECT_SOURCE_DIR}/src/main.cpp)
add_executable(mtxPower ${SOURCES}  ${PROJECT_SOURCE_DIR}/src/matrixPower.cpp)
add_executable(param_study ${SOURCES} ${PROJECT_SOURCE_DIR}/src/param_study.cpp)
add_executable(spmv ${SOURCES}  ${PROJECT_SOURCE_DIR}/src/spmv.cpp)

set(EXT_LIBRARIES "${BLAS_LIBRARIES}")

if(RACE_USE_SPMP)
    set(EXT_LIBRARIES "${PROJECT_BINARY_DIR}/SpMP/src/SpMP/libspmp.a" "${BLAS_LIBRARIES}")
endif()

target_link_libraries(race "${EXT_LIBRARIES}")
target_link_libraries(mtxPower "${EXT_LIBRARIES}")
target_link_libraries(param_study "${EXT_LIBRARIES}")
target_link_libraries(spmv "${EXT_LIBRARIES}")
include_directories(SYSTEM ${PROJECT_BINARY_DIR}/SpMP/src)
#add_dependencies(race SpMP)
#add_dependencies(mtxPower SpMP)
#add_dependencies(param_study SpMP)
#add_dependencies(spmv SpMP)
