cmake_minimum_required(VERSION 3.10)
project(splinelibrary LANGUAGES CXX C)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}")
include_directories(spline_library ${Boost_INCLUDE_DIRS})

add_library(splinelibrary interfaces/spline_library_c.cpp)
target_include_directories(splinelibrary PUBLIC spline_library interfaces)
