cmake_minimum_required(VERSION 3.10)
project(splinelibrary LANGUAGES CXX C)

# use C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(spline_library)

add_library(splinelibrary SHARED interfaces/spline_library_c.cpp)
