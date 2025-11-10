cmake_minimum_required(VERSION 3.27)
project(2dqr)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations")

option(DEBUG "Enable debug mode" OFF)

if (DEBUG)
    add_definitions(-DDEBUG)
endif()

find_package(cxxopts CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)

include_directories(inc)

file(GLOB SRC_COMMON "src/common/*.cpp")
file(GLOB SRC_ALGO "src/algo/*.cpp")
file(GLOB SRC_ABC_INTERFACE "src/abc_interface/*.cpp")

include(ExternalProject)
ExternalProject_Add(abc
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/abc
    BUILD_COMMAND make libabc-pic
    INSTALL_COMMAND ""
)

add_executable(sharp2DQR "src/main.cpp" ${SRC_ABC_INTERFACE} ${SRC_COMMON} ${SRC_ALGO} ${SRC_CUDD_CPP})

add_library(abc-pic STATIC IMPORTED)
set_target_properties(abc-pic PROPERTIES
    IMPORTED_LOCATION ${CMAKE_BINARY_DIR}/abc-prefix/src/abc-build/libabc-pic.a
    INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/abc/src/
)
target_link_libraries(sharp2DQR PRIVATE abc-pic)

target_link_libraries(sharp2DQR PRIVATE cxxopts::cxxopts)
target_link_libraries(sharp2DQR PRIVATE edit)
target_link_libraries(sharp2DQR PRIVATE fmt::fmt)