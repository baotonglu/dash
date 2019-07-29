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


option(USE_PMEM "enable persistent memory support" ON)

set(libs_to_link
    pthread)

if(USE_PMEM MATCHES "ON")
  message(STATUS "persistent memory support enabled, going to build with PMDK")
  add_definitions(-DPMEM)
  list(APPEND libs_to_link pmemobj)
endif()

if(USE_PMEM MATCHES "ON")
  add_executable(test_pmem src/test_pmem.cpp)
  target_link_libraries(test_pmem PRIVATE ${libs_to_link})
else()
  add_executable(test_cuckoo src/test_cuckoo.cpp)
  target_link_libraries(test_cuckoo PRIVATE ${libs_to_link})
endif()
