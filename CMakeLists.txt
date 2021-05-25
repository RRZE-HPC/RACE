cmake_minimum_required (VERSION 2.8)
project (RACE)

set(RACE_VERSION_MAJOR 0)
set(RACE_VERSION_MINOR 3)
set(RACE_VERSION_REVISION 0)
set(RACE_VERSION_BUILD 0)
set(RACE_VERSION_IDENTIFIER)
set(RACE_VERSION "${RACE_VERSION_MAJOR}.${RACE_VERSION_MINOR}.${RACE_VERSION_REVISION}.${RACE_VERSION_BUILD}")

include(ExternalProject)
include(${PROJECT_SOURCE_DIR}/cmake/OptimizeForArchitecture.cmake)

set(RACE_VERBOSITY 1 CACHE STRING "Level of verbosity information. Valid values: 0 (be totally silent), 1 (print warnings, errors and further important information)")
set(RACE_PRINT_ONLYFIRST 1 CACHE BOOL "Print all errors, warnings and info only once for each occurence")

if (NOT BUILD_SHARED_LIBS)
    message(STATUS "No library type specified. Setting it to shared.")
    set(BUILD_SHARED_LIBS ON CACHE BOOL "Set to ON/OFF for shared/static library build.")
endif()


set (SOURCES
    src/graph.cpp
    src/traverse.cpp
    src/levelData.cpp
    src/lb.cpp
    src/error.cpp
    src/zone_tree.cpp
    src/level_recursion.cpp
    src/interface.cpp
    src/functionManager.cpp
    src/machine.cpp
    src/pin.cpp
    src/level_pool.cpp
    src/dummy.cpp
    src/signal.cpp
    #src/simdify.cpp
    src/sell_c_sigmize.cpp
    src/matrixPower.cpp
    src/matrixPowerRecursive.cpp
    )

set (HEADERS
    include/graph.h
    include/graph_SoA.h
    include/graph_AoS.h
    include/traverse.h
    include/traverse_GAP.h
    include/traverse_serial.h
    include/lb.h
    include/error.h
    include/config.h
    include/print.h
    include/type.h
    include/utility.h
    include/interface.h
    include/levelData.h
    include/zone_tree.h
    include/level_recursion.h
    include/functionManager.h
    include/machine.h
    include/pin.h
    include/thpool.h
    include/thpool.tpp
    include/level_pool.h
    include/timing.h
    include/test.h
    include/lock_less.h
    include/fast_barrier.h
    include/spin_cond.h
    include/signal.h
    include/macros.h
    include/simdify.h
    include/sell_c_sigmize.h
    include/matrixPower.h
    include/matrixPowerRecursive.h
    )

set(RACE_PERMUTE_ON_FLY 1 CACHE BOOL "Experimental. Do permutation on fly wthou actually permuting intermediate graph. It speeds up the pre-processing time.")
set(RACE_USE_SOA_GRAPH 1 CACHE BOOL "Experimental. Use SoA graph instead of AoS. It speeds up the pre-processing time.")

set(RACE_USE_GAP 1 CACHE BOOL "Use a modified version of GAP for doing BFS")
if (RACE_USE_GAP)
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
        set(RACE_USE_GAP 1)
        #set(SPMP_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
        set_property(DIRECTORY PROPERTY EP_STEP_TARGETS build)
        ExternalProject_Add(GAP
            PREFIX ${PROJECT_BINARY_DIR}/GAP
            URL https://github.com/sbeamer/gapbs/archive/refs/heads/master.zip
            CONFIGURE_COMMAND ""
            BUILD_COMMAND ""
            BUILD_IN_SOURCE 1
            INSTALL_COMMAND ""
            )

        set(GAP_INCLUDE_DIR "${PROJECT_BINARY_DIR}/GAP/src/GAP/src/")

        set (HEADERS_GAP
            "${PROJECT_BINARY_DIR}/GAP/src/GAP/src/bitmap.h"
            "${PROJECT_BINARY_DIR}/GAP/src/GAP/src/platform_atomics.h"
            "${PROJECT_BINARY_DIR}/GAP/src/GAP/src/sliding_queue.h"
            )
    else()
        message(FATAL_ERROR "GAP BFS not used. Pre-processing will be slow")
    endif()
endif()



set(HWLOC_INCLUDE_DIR HWLOC_INCLUDE_DIR-NOTFOUND CACHE PATH "Location of hwloc.h")
set(LIBHWLOC LIBHWLOC-NOTFOUND)

if (PKGCONFIG_FOUND AND NOT HWLOC_INCLUDE_DIR)
    message(STATUS "Try to find hwloc via pkg-config")
    find_package(HWLOC)
    find_library(LIBHWLOC hwloc)
endif()

# find hwloc.h in default header locations
if (NOT HWLOC_INCLUDE_DIR)
    find_path(HWLOC_INCLUDE_DIR hwloc.h)
endif()

#set(GITHUB_GIT_URL "github.com" CACHE STRING "The github GIT URL. This has to be changed in certain situations if located behind a firewall.")
set(BUILTIN_HWLOC false)

# hwloc not found - download and build it
if (NOT HWLOC_INCLUDE_DIR)
    message(STATUS "Could not find hwloc.h neither via pkg-config nor in the defined HWLOC_INCLUDE_DIR! Trying to build hwloc...")
    ExternalProject_Add(hwloc
        PREFIX ${PROJECT_BINARY_DIR}/hwloc
        URL https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.6.tar.gz
        SOURCE_DIR ${PROJECT_BINARY_DIR}/hwloc/src/hwloc
        #CONFIGURE_COMMAND ./configure --enable-static --disable-shared --prefix=${PROJECT_BINARY_DIR}/hwloc-inst
        CONFIGURE_COMMAND CXXFLAGS=-fPIC CFLAGS=-fPIC ./configure --prefix=${PROJECT_BINARY_DIR}/hwloc-inst
        BUILD_COMMAND make
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND make install
        )
    set(BUILTIN_HWLOC 1)
    set(HWLOC_INCLUDE_DIR "${PROJECT_BINARY_DIR}/hwloc-inst/include")
    set(LIBHWLOC "${PROJECT_BINARY_DIR}/hwloc-inst/lib/libhwloc.so")
    #set(LIBHWLOC "${PROJECT_BINARY_DIR}/hwloc-inst/lib64/libhwloc.so")
endif()

if (NOT LIBHWLOC)
    get_filename_component(HWLOC_LIBDIR ${HWLOC_INCLUDE_DIR}/../lib REALPATH)
    find_library(LIBHWLOC hwloc HINTS ${HWLOC_LIBDIR})
endif()

set(HWLOC_LIBRARIES ${LIBHWLOC})

set(RACE_THREAD_METHOD OpenMP CACHE STRING "Which threading method should be used for the kernel: options are OpenMP or PThread")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fPIC -std=c99 -Wall -Wextra -Werror -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fPIC -std=c++0x -Wall -Wextra -Werror -pedantic")
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_INCLUDE_SYSTEM_FLAG_C "-isystem ") # Not present in Cmake 2.8.11.2
    set(CMAKE_INCLUDE_SYSTEM_FLAG_CXX "-isystem ") # Not present in Cmake 2.8.11.2
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wunused-variable -save-temps -fsource-asm -w3 -wd981")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunused-variable -save-temps -fsource-asm -w3 -wd981 -wd383")
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fno-alias")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-alias")
    endif()
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-pragmas")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") 
        message(STATUS "NOTE: by default, we only set the -save-temps compiler flag with Intel compilers right now. You may want to use CC=\"ccache gcc\" etc to speed up the compiliation.")
        if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -fstrict-aliasing")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -fstrict-aliasing")
        endif()
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT "${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "4.8"))
            set(SANITIZER "" CACHE STRING "What should be sanitized by Clang/GCC for a build? (address, thread, ..., undefined or leave empty)")
            if (NOT "${SANITIZER}" STREQUAL "") 
                string(TOLOWER "${SANITIZER}" SANITIZER)
                set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -fsanitize=${SANITIZER} -fno-omit-frame-pointer")
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -fsanitize=${SANITIZER} -fno-omit-frame-pointer")
                if (GHOST_HAVE_CUDA)
                    message(WARNING "Sanitizers are probably not working together with CUDA!")
                endif()
            endif()
        endif()
    endif()
endif()

find_package(OpenMP)

if (OPENMP_FOUND)
    set(RACE_USE_OPENMP 1 CACHE BOOL "Build with OpenMP")
    set(RACE_HAVE_OPENMP 1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else ()
    set(RACE_USE_OPENMP 0 CACHE BOOL "Build with OpenMP")
    set(RACE_HAVE_OPENMP 0)
endif()

if(NOT "${RACE_THREAD_METHOD}" STREQUAL "PThread")
    set(RACE_KERNEL_THREAD_OMP 1)
else ()
    set(RACE_KERNEL_THREAD_OMP 0)
endif()

configure_file (
    "${PROJECT_SOURCE_DIR}/include/config.h.in"
    "${PROJECT_SOURCE_DIR}/include/config.h"
)
configure_file(${PROJECT_SOURCE_DIR}/race-config.cmake.in ${PROJECT_BINARY_DIR}/race-config.cmake @ONLY)

include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(RACE ${SOURCES})

#get_filename_component(HWLOC_LIB_PATH ${LIBHWLOC} DIRECTORY)
set(RACE_LIB_PATH ${CMAKE_INSTALL_PREFIX}/lib/RACE)

#if (${BUILD_SHARED_LIBS})
#    file(WRITE ${PROJECT_BINARY_DIR}/RACE_LFlags.sh "#!/bin/bash\n echo \"-L ${RACE_LIB_PATH} -Wl,-rpath ${RACE_LIB_PATH} -L ${HWLOC_LIB_PATH} -Wl,-rpath ${HWLOC_LIB_PATH} -l hwloc -l RACE\"")
#else ()
#    file(WRITE ${PROJECT_BINARY_DIR}/RACE_LFlags.sh "#!/bin/bash\n echo \"${RACE_LIB_PATH}/libRACE.a -L ${HWLOC_LIB_PATH} -Wl,-rpath ${HWLOC_LIB_PATH} -l hwloc\"")
#endif()

if(BUILTIN_HWLOC)
    target_link_libraries(RACE "${LIBHWLOC}")
    include_directories(SYSTEM "${HWLOC_INCLUDE_DIR}")
    add_dependencies(RACE hwloc)
endif()

install(FILES ${HEADERS} DESTINATION include/RACE/)
install(FILES ${HEADERS_GAP} DESTINATION include/RACE/)

#install(FILES ${PROJECT_BINARY_DIR}/RACE_LFlags.sh DESTINATION lib/RACE PERMISSIONS OWNER_READ OWNER_EXECUTE)
install(FILES ${PROJECT_BINARY_DIR}/race-config.cmake DESTINATION lib/RACE)

if (${BUILD_SHARED_LIBS})
    install(TARGETS RACE EXPORT RACE-targets LIBRARY DESTINATION lib/RACE)
else ()
    install(TARGETS RACE EXPORT RACE-targets ARCHIVE DESTINATION lib/RACE)
endif()


## BUILD perf_tests ##

#AutodetectHostArchitecture()
OptimizeForArchitecture()

#Set Vector length
if(USE_AVX512F)
    set(VECTOR_LENGTH 8)
elseif(USE_AVX OR USE_AVX2)
    set(VECTOR_LENGTH 4)
elseif(USE_SSE2 OR USE_SSE3 OR SSE4_1 OR SSE4_2 OR SSE_4a)
    set(VECTOR_LENGTH 2)
else()
    set(VECTOR_LENGTH 1)
endif()


set(PERF_TEST_USE_SPMP 1 CACHE BOOL "Use RCM by SpMP for matrix preprocessing")
if (PERF_TEST_USE_SPMP)
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


configure_file (
    "${PROJECT_SOURCE_DIR}/perf_test/include/config_eg.h.in"
    "${PROJECT_BINARY_DIR}/perf_test/include/config_eg.h"
    )

include_directories(${PROJECT_SOURCE_DIR}/perf_test/include)
include_directories(${PROJECT_BINARY_DIR}/perf_test/include)

file(GLOB PERF_SOURCES ${PROJECT_SOURCE_DIR}/perf_test/src/helpers/*.cpp)

add_executable(race ${PERF_SOURCES}  ${PROJECT_SOURCE_DIR}/perf_test/src/race.cpp)
add_executable(mtxPower ${PERF_SOURCES}  ${PROJECT_SOURCE_DIR}/perf_test/src/matrixPower.cpp)

set_target_properties(race PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/perf_test")
set_target_properties(mtxPower PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/perf_test")

if(PERF_TEST_USE_SPMP)
    set(EXT_LIBRARIES "${PROJECT_BINARY_DIR}/SpMP/src/SpMP/libspmp.a")
endif()

if(RACE_USE_GAP)
    include_directories(SYSTEM "${GAP_INCLUDE_DIR}")
    add_dependencies(RACE GAP)
endif()


if(BUILTIN_HWLOC)
    target_link_libraries(race "${EXT_LIBRARIES}" RACE)
    target_link_libraries(mtxPower "${EXT_LIBRARIES}" RACE )
    add_dependencies(race hwloc)
    add_dependencies(mtxPower hwloc)
else()
    target_link_libraries(race "${EXT_LIBRARIES}" RACE hwloc)
    target_link_libraries(mtxPower "${EXT_LIBRARIES}" RACE hwloc)
endif()

add_dependencies(race RACE)
add_dependencies(mtxPower RACE)
if(PERF_TEST_USE_SPMP)
    include_directories(SYSTEM ${PROJECT_BINARY_DIR}/SpMP/src)
    add_dependencies(race SpMP)
    add_dependencies(mtxPower SpMP)
endif()

file(COPY ${PROJECT_SOURCE_DIR}/perf_test/perf_script.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/perf_test/.)
file(COPY ${PROJECT_SOURCE_DIR}/perf_test/check-state.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/perf_test/.)
file(COPY ${PROJECT_SOURCE_DIR}/perf_test/machine-state.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/perf_test/.)
file(COPY ${PROJECT_SOURCE_DIR}/perf_test/config.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/perf_test/.)
file(COPY ${PROJECT_SOURCE_DIR}/perf_test/config.txt DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/perf_test/.)
file(COPY ${PROJECT_SOURCE_DIR}/perf_test/plot_generator DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/perf_test/.)
