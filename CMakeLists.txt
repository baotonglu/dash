if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif ()

project(LinearHashing)
cmake_minimum_required(VERSION 3.11)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-long-long -pedantic -fPIC")
set(CMAKE_ENABLE_COMPILE_COMMANDS "ON")

if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
  message("-- Release mode, all optimizations enabled")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Ofast -march=native")
endif ()


add_executable(test_cuckoo src/test_cuckoo.cpp)
target_include_directories(test_cuckoo PRIVATE src util)
target_link_libraries(test_cuckoo pthread)
