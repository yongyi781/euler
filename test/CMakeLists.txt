cmake_minimum_required(VERSION 3.30)
project(euler VERSION 0.1.0 LANGUAGES CXX)
include(CTest)

find_package(primecount REQUIRED)
find_package(TBB)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_SCAN_FOR_MODULES OFF)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(-march=native)
endif()

include_directories(../include)
link_libraries(primecount::primecount primesieve::primesieve gmp)

if (TBB_FOUND)
    link_libraries(TBB::tbb)
endif()

file(GLOB sources *.cpp)

foreach(source_file ${sources})
    cmake_path(GET source_file STEM target)
    message("Adding test ${target}")
    add_executable(${target} ${source_file})
    add_test(NAME ${target} COMMAND ${target})
endforeach()
