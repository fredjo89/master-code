project(mpitry)

cmake_minimum_required(VERSION 2.6)

# We use the C++ language
enable_language(CXX)

# Include MPI
find_package(MPI)
add_definitions(-DHAVE_MPI=1)
set(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
set(DEPLIBS_CXX ${DEPLIBS_CXX} ${MPI_CXX_LIBRARIES})

# Send the include directories to the compiler
include_directories(${INCLUDES})

#SET(CMAKE_BUILD_TYPE distribution)
#SET(CMAKE_CXX_FLAGS_DISTRIBUTION "-O1")


set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

# Add libaries
add_library(read readFromFile.cpp)
add_library(dataDefs data_defs.cpp)
add_library(utils utilities.cpp)
add_library(distr distribution.cpp)
add_library(jacobi jacobi.cpp)

# Add executable
add_executable(main main.cpp)

# Link to libraries
target_link_libraries(main ${DEPLIBS_CXX})
target_link_libraries (main read dataDefs distr utils jacobi -lmetis )
# We want to link this as a C++ program
set_target_properties(main PROPERTIES LINKER_LANGUAGE CXX)
