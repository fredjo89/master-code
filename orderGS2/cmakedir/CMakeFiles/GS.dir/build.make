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
CMAKE_SOURCE_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/cmakedir

# Include any dependencies generated for this target.
include CMakeFiles/GS.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/GS.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/GS.dir/flags.make

CMakeFiles/GS.dir/orderGS_Solver.cpp.o: CMakeFiles/GS.dir/flags.make
CMakeFiles/GS.dir/orderGS_Solver.cpp.o: ../orderGS_Solver.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/cmakedir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/GS.dir/orderGS_Solver.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/GS.dir/orderGS_Solver.cpp.o -c /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/orderGS_Solver.cpp

CMakeFiles/GS.dir/orderGS_Solver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/GS.dir/orderGS_Solver.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/orderGS_Solver.cpp > CMakeFiles/GS.dir/orderGS_Solver.cpp.i

CMakeFiles/GS.dir/orderGS_Solver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/GS.dir/orderGS_Solver.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/orderGS_Solver.cpp -o CMakeFiles/GS.dir/orderGS_Solver.cpp.s

CMakeFiles/GS.dir/orderGS_Solver.cpp.o.requires:
.PHONY : CMakeFiles/GS.dir/orderGS_Solver.cpp.o.requires

CMakeFiles/GS.dir/orderGS_Solver.cpp.o.provides: CMakeFiles/GS.dir/orderGS_Solver.cpp.o.requires
	$(MAKE) -f CMakeFiles/GS.dir/build.make CMakeFiles/GS.dir/orderGS_Solver.cpp.o.provides.build
.PHONY : CMakeFiles/GS.dir/orderGS_Solver.cpp.o.provides

CMakeFiles/GS.dir/orderGS_Solver.cpp.o.provides.build: CMakeFiles/GS.dir/orderGS_Solver.cpp.o

# Object files for target GS
GS_OBJECTS = \
"CMakeFiles/GS.dir/orderGS_Solver.cpp.o"

# External object files for target GS
GS_EXTERNAL_OBJECTS =

libGS.a: CMakeFiles/GS.dir/orderGS_Solver.cpp.o
libGS.a: CMakeFiles/GS.dir/build.make
libGS.a: CMakeFiles/GS.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libGS.a"
	$(CMAKE_COMMAND) -P CMakeFiles/GS.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/GS.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/GS.dir/build: libGS.a
.PHONY : CMakeFiles/GS.dir/build

CMakeFiles/GS.dir/requires: CMakeFiles/GS.dir/orderGS_Solver.cpp.o.requires
.PHONY : CMakeFiles/GS.dir/requires

CMakeFiles/GS.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/GS.dir/cmake_clean.cmake
.PHONY : CMakeFiles/GS.dir/clean

CMakeFiles/GS.dir/depend:
	cd /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/cmakedir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2 /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2 /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/orderGS2/cmakedir/CMakeFiles/GS.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/GS.dir/depend

