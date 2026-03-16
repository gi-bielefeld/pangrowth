cmake_minimum_required(VERSION 3.5)
project(pangrowth)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

#set(CMAKE_CXX_FLAGS "-Wall -Wextra -O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

file(GLOB HEAD_C "src/*.h")
file(GLOB SRC_C "src/*.c")
#file(GLOB HEAD_CXX "src/*.hpp")
file(GLOB SRC_CXX "src/*.cpp")

add_executable(pangrowth ${HEAD_C} ${SRC_C} ${HEAD_CXX} ${SRC_CXX})

find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
  target_link_libraries(pangrowth OpenMP::OpenMP_CXX)
endif()

find_package(Rmath 1.0.0 REQUIRED)
target_link_libraries(pangrowth Rmath::Rmath z m pthread)
#target_link_libraries(pangrowth z m pthread)
