cmake_minimum_required(VERSION 3.1)
project(assignment3)

if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "No build type selected, default to Release")
	set(CMAKE_BUILD_TYPE "Release")
endif()

# Project sources
add_executable(${PROJECT_NAME}
	src/main.cpp
	src/utils.h
)

# Include Eigen for linear algebra, stb and gif-h to export images
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../ext/eigen" "${CMAKE_CURRENT_SOURCE_DIR}/../ext/stb" "${CMAKE_CURRENT_SOURCE_DIR}/../ext/gif-h")

# Use C++11 version of the standard
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON)

# Place the output binary at the root of the build folder
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
