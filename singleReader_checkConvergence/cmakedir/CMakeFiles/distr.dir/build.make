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
CMAKE_SOURCE_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/cmakedir

# Include any dependencies generated for this target.
include CMakeFiles/distr.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/distr.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/distr.dir/flags.make

CMakeFiles/distr.dir/distribution.cpp.o: CMakeFiles/distr.dir/flags.make
CMakeFiles/distr.dir/distribution.cpp.o: ../distribution.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/cmakedir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/distr.dir/distribution.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/distr.dir/distribution.cpp.o -c /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/distribution.cpp

CMakeFiles/distr.dir/distribution.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/distr.dir/distribution.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/distribution.cpp > CMakeFiles/distr.dir/distribution.cpp.i

CMakeFiles/distr.dir/distribution.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/distr.dir/distribution.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/distribution.cpp -o CMakeFiles/distr.dir/distribution.cpp.s

CMakeFiles/distr.dir/distribution.cpp.o.requires:
.PHONY : CMakeFiles/distr.dir/distribution.cpp.o.requires

CMakeFiles/distr.dir/distribution.cpp.o.provides: CMakeFiles/distr.dir/distribution.cpp.o.requires
	$(MAKE) -f CMakeFiles/distr.dir/build.make CMakeFiles/distr.dir/distribution.cpp.o.provides.build
.PHONY : CMakeFiles/distr.dir/distribution.cpp.o.provides

CMakeFiles/distr.dir/distribution.cpp.o.provides.build: CMakeFiles/distr.dir/distribution.cpp.o

# Object files for target distr
distr_OBJECTS = \
"CMakeFiles/distr.dir/distribution.cpp.o"

# External object files for target distr
distr_EXTERNAL_OBJECTS =

libdistr.a: CMakeFiles/distr.dir/distribution.cpp.o
libdistr.a: CMakeFiles/distr.dir/build.make
libdistr.a: CMakeFiles/distr.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libdistr.a"
	$(CMAKE_COMMAND) -P CMakeFiles/distr.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/distr.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/distr.dir/build: libdistr.a
.PHONY : CMakeFiles/distr.dir/build

CMakeFiles/distr.dir/requires: CMakeFiles/distr.dir/distribution.cpp.o.requires
.PHONY : CMakeFiles/distr.dir/requires

CMakeFiles/distr.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/distr.dir/cmake_clean.cmake
.PHONY : CMakeFiles/distr.dir/clean

CMakeFiles/distr.dir/depend:
	cd /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/cmakedir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/singleReader_checkConvergence/cmakedir/CMakeFiles/distr.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/distr.dir/depend
