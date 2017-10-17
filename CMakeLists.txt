cmake_minimum_required (VERSION 2.8)
project (RACE)

set(RACE_VERSION_MAJOR 0)
set(RACE_VERSION_MINOR 0)
set(RACE_VERSION_REVISION 0)
set(RACE_VERSION_IDENTIFIER)
set(RACE_VERSION "${RACE_VERSION_MAJOR}.${RACE_VERSION_MINOR}.${RACE_VERSION_REVISION}")

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
    )

set (HEADERS
    include/graph.h
    include/traverse.h
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
    )


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

# hwloc not found - download and build it
if (NOT HWLOC_INCLUDE_DIR)
    message(STATUS "Could not find hwloc.h neither via pkg-config nor in the defined HWLOC_INCLUDE_DIR! Trying to build hwloc...")
    ExternalProject_Add(hwloc
        PREFIX ${PROJECT_BINARY_DIR}/hwloc
        URL https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.6.tar.gz
        SOURCE_DIR ${PROJECT_BINARY_DIR}/hwloc/src/hwloc
        CONFIGURE_COMMAND CXXFLAGS=-fPIC CFLAGS=-fPIC ./configure --enable-static --disable-shared --prefix=${PROJECT_BINARY_DIR}/hwloc-inst
        BUILD_COMMAND make
        BUILD_IN_SOURCE 1
        INSTALL_COMMAND make install
        )
    set(BUILTIN_HWLOC 1)
    set(HWLOC_INCLUDE_DIR "${PROJECT_BINARY_DIR}/hwloc-inst/include")
    set(LIBHWLOC "${PROJECT_BINARY_DIR}/hwloc-inst/lib/libhwloc.a")
endif()

if (NOT LIBHWLOC)
    get_filename_component(HWLOC_LIBDIR ${HWLOC_INCLUDE_DIR}/../lib REALPATH)
    find_library(LIBHWLOC hwloc HINTS ${HWLOC_LIBDIR})
endif()

set(HWLOC_LIBRARIES ${LIBHWLOC})

set(RACE_THREAD_METHOD PThread CACHE STRING "Which threading method should be used for the kernel: options are OpenMP or PThread")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fPIC -std=c99 -Wall -Wextra -Werror -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -std=c++0x -Wall -Wextra -Werror -pedantic")
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
        message(STATUS "NOTE: by default, we only set the -save-temps compiler flag with Intel compilers right now. You may want to use CC=\"ccache gcc\" etc to speed up the compiliation, see issue #269.")
        if ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Ofast -fstrict-aliasing")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -fstrict-aliasing")
        endif()
        if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT "${CMAKE_CXX_COMPILER_VERSION}" VERSION_LESS "4.8"))
            set(SANITIZER "" CACHE STRING "What should be sanitized by Clang/GCC for a build? (address, thread, memory, undefined or leave empty)")
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


include_directories(${PROJECT_SOURCE_DIR}/include)
add_library(RACE ${SOURCES})

target_link_libraries(RACE ${LIBHWLOC})
include_directories(SYSTEM ${HWLOC_INCLUDE_DIR})
add_dependencies(RACE hwloc)

install(FILES ${HEADERS} DESTINATION include/RACE/)
if (${BUILD_SHARED_LIBS})
    install(TARGETS RACE EXPORT RACE-targets LIBRARY DESTINATION lib/RACE)
else ()
    install(TARGETS RACE EXPORT RACE-targets ARCHIVE DESTINATION lib/RACE)
endif()

