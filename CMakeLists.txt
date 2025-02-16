cmake_minimum_required(VERSION 3.10)

project(p-sort-index VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(ENABLE_FEATURE_X "Enable Feature X" OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(p-sort-index ${SRC_FILES})

target_compile_options(p-sort-index PRIVATE -O3) # O3 and NDEBUG

find_package(Boost 1.70 REQUIRED COMPONENTS date_time filesystem iostreams graph program_options)

link_directories(~/lib)

# Link libraries (equivalent to -lsdsl -ldivsufsort -ldivsufsort64)
target_link_libraries(p-sort-index PRIVATE sdsl divsufsort divsufsort64)

find_package(Threads REQUIRED)

# Link pthread
target_link_libraries(p-sort-index PRIVATE Threads::Threads  Boost::program_options 
Boost::iostreams 
Boost::graph
stdc++)