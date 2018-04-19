cmake_minimum_required(VERSION 2.8)

# project name
project(RACE_example CXX C)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_PREFIX_PATH ${RACE_DIR})
# find RACE
find_package(RACE REQUIRED CONFIG)
message(STATUS "Found RACE")

include_directories(${RACE_INCLUDE_DIR})
link_libraries(${RACE_LIBRARIES})
set(CXX_STANDARD 11)
Set(CXX_STANDARD_REQUIRED On)
include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${RACE_INCLUDE_DIRS})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${RACE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} ${RACE_CXX_FLAGS}")

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -xHOST -save-temps -fsource-asm")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xHOST -save-temps -fsource-asm")
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

file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)
add_executable(race ${SOURCES})
