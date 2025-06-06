cmake_minimum_required(VERSION 3.30)
project(euler VERSION 0.1.0 LANGUAGES CXX)
find_package(TBB)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

include_directories(../include SYSTEM C:/Tools/primecount/include)
if(CMAKE_CXX_PLATFORM_ID STREQUAL "Windows")
    message("Adding Windows-specific headers")
    include_directories(SYSTEM C:/Tools/boost_1_84_0 C:/Tools/gmp/include)
    link_directories(SYSTEM C:/Tools/primecount-msvc/lib)
else()
    link_directories(SYSTEM C:/Tools/primecount/lib)
    link_libraries(stdc++exp TBB::tbb gmp)
endif()
link_libraries(primecount primesieve)

include(CTest)
enable_testing()

file(GLOB sources *.cpp)

foreach(source_file ${sources})
    cmake_path(GET source_file STEM target)
    message("Adding test ${target}")
    add_executable(${target} ${source_file})
    add_test(NAME ${target} COMMAND ${target})
endforeach()
