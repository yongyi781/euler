cmake_minimum_required(VERSION 3.30)
project(euler VERSION 0.1.0 LANGUAGES CXX)
find_package(TBB)

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

include_directories(../include SYSTEM C:/Tools/primecount/include)
if(CMAKE_CXX_PLATFORM_ID STREQUAL "Windows")
    message("Adding Windows-specific headers")
    include_directories(SYSTEM C:/Tools/boost_1_84_0 C:/Tools/gmp/include)
else()
    link_libraries(stdc++exp TBB::tbb)
endif()
link_directories(SYSTEM C:/Tools/primecount/lib)
link_libraries(primecount primesieve gmp)

include(CTest)
enable_testing()

set(targets
    basic
    crt
    find_recurrence)

foreach(target ${targets})
    message("Adding target: ${target}")
    add_executable(${target} "${target}.cpp")
    add_test(NAME ${target} COMMAND ${target})
endforeach()
