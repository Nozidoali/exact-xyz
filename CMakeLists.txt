# set required cmake version
cmake_minimum_required(VERSION 3.20)

project(
    xyz
    LANGUAGES CXX
    VERSION 0.1.0
    DESCRIPTION "quantum computing library"
)

option(BUILD_XYZ_BINDINGS "Build the XYZ Python bindings" OFF)
option(BUILD_XYZ_TESTS "Build the XYZ tests" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

set(XYZ_SOURCE_DIR ${CMAKE_SOURCE_DIR})

set(PACKAGE_NAME cxyz)


find_package(fmt REQUIRED)
get_target_property(FMT_INCLUDE_DIRS fmt::fmt INTERFACE_INCLUDE_DIRECTORIES)
include_directories(${FMT_INCLUDE_DIRS})
message(STATUS "FMT_INCLUDE_DIRS: ${FMT_INCLUDE_DIRS}")

include_directories(${XYZ_SOURCE_DIR}/include)

add_subdirectory(${XYZ_SOURCE_DIR}/lib)
add_subdirectory(${XYZ_SOURCE_DIR}/tools)

target_link_libraries(${PACKAGE_NAME} PUBLIC fmt::fmt)

if (BUILD_XYZ_BINDINGS)
    add_subdirectory(${XYZ_SOURCE_DIR}/bindings)
endif()

if (BUILD_XYZ_TESTS)
    add_subdirectory(${XYZ_SOURCE_DIR}/tests)
endif()