project(mpitry)

cmake_minimum_required(VERSION 2.6)

# We use the C and C++ language
enable_language(C)
enable_language(CXX)

# Include MPI
find_package(MPI)
add_definitions(-DHAVE_MPI=1)
set(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
set(DEPLIBS_CXX ${DEPLIBS_CXX} ${MPI_CXX_LIBRARIES})
set(DEPLIBS_F ${DEPLIBS_F} ${MPI_Fortran_LIBRARIES})
# Send the include directories to the compiler
include_directories(${INCLUDES})


set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE Debug)

# Add basis_solver.cpp to library


# Add executable
add_executable(main main.cpp)
add_executable(main2 main2.cpp)

# Link to libraries
target_link_libraries(main ${DEPLIBS_CXX})
target_link_libraries (main   )
target_link_libraries(main2 ${DEPLIBS_CXX})
# We want to link this as a C program
set_target_properties(main PROPERTIES LINKER_LANGUAGE CXX)
