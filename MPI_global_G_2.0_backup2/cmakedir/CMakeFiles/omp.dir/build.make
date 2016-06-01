# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/cmakedir

# Include any dependencies generated for this target.
include CMakeFiles/omp.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/omp.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/omp.dir/flags.make

CMakeFiles/omp.dir/ompSolver.cpp.o: CMakeFiles/omp.dir/flags.make
CMakeFiles/omp.dir/ompSolver.cpp.o: ../ompSolver.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/cmakedir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/omp.dir/ompSolver.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/omp.dir/ompSolver.cpp.o -c /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/ompSolver.cpp

CMakeFiles/omp.dir/ompSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/omp.dir/ompSolver.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/ompSolver.cpp > CMakeFiles/omp.dir/ompSolver.cpp.i

CMakeFiles/omp.dir/ompSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/omp.dir/ompSolver.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/ompSolver.cpp -o CMakeFiles/omp.dir/ompSolver.cpp.s

CMakeFiles/omp.dir/ompSolver.cpp.o.requires:
.PHONY : CMakeFiles/omp.dir/ompSolver.cpp.o.requires

CMakeFiles/omp.dir/ompSolver.cpp.o.provides: CMakeFiles/omp.dir/ompSolver.cpp.o.requires
	$(MAKE) -f CMakeFiles/omp.dir/build.make CMakeFiles/omp.dir/ompSolver.cpp.o.provides.build
.PHONY : CMakeFiles/omp.dir/ompSolver.cpp.o.provides

CMakeFiles/omp.dir/ompSolver.cpp.o.provides.build: CMakeFiles/omp.dir/ompSolver.cpp.o

# Object files for target omp
omp_OBJECTS = \
"CMakeFiles/omp.dir/ompSolver.cpp.o"

# External object files for target omp
omp_EXTERNAL_OBJECTS =

libomp.a: CMakeFiles/omp.dir/ompSolver.cpp.o
libomp.a: CMakeFiles/omp.dir/build.make
libomp.a: CMakeFiles/omp.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libomp.a"
	$(CMAKE_COMMAND) -P CMakeFiles/omp.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/omp.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/omp.dir/build: libomp.a
.PHONY : CMakeFiles/omp.dir/build

CMakeFiles/omp.dir/requires: CMakeFiles/omp.dir/ompSolver.cpp.o.requires
.PHONY : CMakeFiles/omp.dir/requires

CMakeFiles/omp.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/omp.dir/cmake_clean.cmake
.PHONY : CMakeFiles/omp.dir/clean

CMakeFiles/omp.dir/depend:
	cd /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/cmakedir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0 /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0 /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/MPI_global_G_2.0/cmakedir/CMakeFiles/omp.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/omp.dir/depend

