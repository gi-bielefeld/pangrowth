cmake_minimum_required(VERSION 3.10)
project(pangrowth)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#set(CMAKE_CXX_FLAGS "-Wall -Wextra -O2")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

file(GLOB HEAD_C "src/*.h")
file(GLOB SRC_C "src/*.c")
#file(GLOB HEAD_CXX "src/*.hpp")
file(GLOB SRC_CXX "src/*.cpp")

option(PANGROWTH_WITH_GGCAT "Enable ggcat C++ API integration" OFF)
set(GGCAT_CPP_API_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/ggcat-cpp-api" CACHE PATH "Path to ggcat-cpp-api wrapper")
set(GGCAT_CPP_BINDINGS_LIB "${GGCAT_CPP_API_DIR}/lib/libggcat_cpp_bindings.a" CACHE FILEPATH "Path to libggcat_cpp_bindings.a")
set(GGCAT_CXX_INTEROP_LIB "${GGCAT_CPP_API_DIR}/lib/libggcat_cxx_interop.a" CACHE FILEPATH "Path to libggcat_cxx_interop.a")

if(PANGROWTH_WITH_GGCAT)
  if(NOT EXISTS "${GGCAT_CPP_BINDINGS_LIB}")
    message(FATAL_ERROR "PANGROWTH_WITH_GGCAT requires GGCAT_CPP_BINDINGS_LIB=/path/to/libggcat_cpp_bindings.a")
  endif()
  if(NOT EXISTS "${GGCAT_CXX_INTEROP_LIB}")
    message(FATAL_ERROR "PANGROWTH_WITH_GGCAT requires GGCAT_CXX_INTEROP_LIB=/path/to/libggcat_cxx_interop.a")
  endif()
  list(APPEND SRC_CXX "${GGCAT_CPP_API_DIR}/src/ggcat.cc")
endif()

add_executable(pangrowth ${HEAD_C} ${SRC_C} ${HEAD_CXX} ${SRC_CXX})

if(PANGROWTH_WITH_GGCAT)
  target_compile_definitions(pangrowth PRIVATE PANGROWTH_WITH_GGCAT=1)
  target_include_directories(pangrowth PRIVATE "${GGCAT_CPP_API_DIR}/include")
  target_link_libraries(pangrowth ${GGCAT_CPP_BINDINGS_LIB} ${GGCAT_CXX_INTEROP_LIB} dl)
endif()

find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
  target_link_libraries(pangrowth OpenMP::OpenMP_CXX)
endif()

find_package(Rmath 1.0.0 REQUIRED)
target_link_libraries(pangrowth Rmath::Rmath z m pthread)
#target_link_libraries(pangrowth z m pthread)
