cmake_minimum_required(VERSION 3.1)
project(Assignment_6)


################################################################################

# Directory to external libraries used in the project
set(THIRD_PARTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../ext/)

# Project sources
add_library(${PROJECT_NAME}
	src/SDLViewer.cpp
	src/SDLViewer.h
)

# SDL
target_compile_definitions(${PROJECT_NAME} PUBLIC SDL_MAIN_HANDLED)
add_subdirectory(${THIRD_PARTY_DIR}/SDL ${CMAKE_BINARY_DIR}/SDL-build)
target_link_libraries(${PROJECT_NAME} PUBLIC SDL2-static)

# if(WIN32)
#     target_link_libraries(${PROJECT_NAME} PUBLIC user32 gdi32 winmm imm32 ole32 oleaut32 version uuid advapi32 setupapi shell32 vcruntime)
# endif()

################################################################################

# Use C++11 version of the standard
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED YES)

target_compile_features(${PROJECT_NAME}
    PUBLIC
        cxx_nullptr
        cxx_lambdas
)

# Place the output binary at the root of the build folder
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# Include Eigen for linear algebra, stb and gif-h for exporting images
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC "${THIRD_PARTY_DIR}/eigen" "${THIRD_PARTY_DIR}/stb" "${THIRD_PARTY_DIR}/gif-h")


################################################################################
################################################################################

add_executable(RasterViewer src/raster.cpp src/RasterViewer.cpp)
target_link_libraries(RasterViewer PUBLIC ${PROJECT_NAME})

# Folder where data files are stored (meshes & stuff)
set(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/")
target_compile_definitions(RasterViewer PUBLIC -DDATA_DIR=\"${DATA_DIR}\")
