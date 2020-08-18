cmake_minimum_required(VERSION 3.1)
project(assignment4)

# Project sources
add_executable(${PROJECT_NAME}
	src/main.cpp
	src/utils.h
)

# Include Eigen for linear algebra, stb and gif-h to export images, json to parse the json files
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../ext/eigen" "${CMAKE_CURRENT_SOURCE_DIR}/../ext/stb" "${CMAKE_CURRENT_SOURCE_DIR}/../ext/gif-h" "${CMAKE_CURRENT_SOURCE_DIR}/../ext/json")

# Use C++11 version of the standard
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON)

# Place the output binary at the root of the build folder
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# Folder where data files are stored (meshes & stuff)
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/")
target_compile_definitions(${PROJECT_NAME} PUBLIC -DDATA_DIR=\"${DATA_DIR}\")
