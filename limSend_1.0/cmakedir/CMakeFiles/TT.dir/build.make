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
CMAKE_SOURCE_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/cmakedir

# Include any dependencies generated for this target.
include CMakeFiles/TT.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/TT.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TT.dir/flags.make

CMakeFiles/TT.dir/TypeTwo.cpp.o: CMakeFiles/TT.dir/flags.make
CMakeFiles/TT.dir/TypeTwo.cpp.o: ../TypeTwo.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/cmakedir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/TT.dir/TypeTwo.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/TT.dir/TypeTwo.cpp.o -c /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/TypeTwo.cpp

CMakeFiles/TT.dir/TypeTwo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TT.dir/TypeTwo.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/TypeTwo.cpp > CMakeFiles/TT.dir/TypeTwo.cpp.i

CMakeFiles/TT.dir/TypeTwo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TT.dir/TypeTwo.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/TypeTwo.cpp -o CMakeFiles/TT.dir/TypeTwo.cpp.s

CMakeFiles/TT.dir/TypeTwo.cpp.o.requires:
.PHONY : CMakeFiles/TT.dir/TypeTwo.cpp.o.requires

CMakeFiles/TT.dir/TypeTwo.cpp.o.provides: CMakeFiles/TT.dir/TypeTwo.cpp.o.requires
	$(MAKE) -f CMakeFiles/TT.dir/build.make CMakeFiles/TT.dir/TypeTwo.cpp.o.provides.build
.PHONY : CMakeFiles/TT.dir/TypeTwo.cpp.o.provides

CMakeFiles/TT.dir/TypeTwo.cpp.o.provides.build: CMakeFiles/TT.dir/TypeTwo.cpp.o

# Object files for target TT
TT_OBJECTS = \
"CMakeFiles/TT.dir/TypeTwo.cpp.o"

# External object files for target TT
TT_EXTERNAL_OBJECTS =

libTT.a: CMakeFiles/TT.dir/TypeTwo.cpp.o
libTT.a: CMakeFiles/TT.dir/build.make
libTT.a: CMakeFiles/TT.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libTT.a"
	$(CMAKE_COMMAND) -P CMakeFiles/TT.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/TT.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TT.dir/build: libTT.a
.PHONY : CMakeFiles/TT.dir/build

CMakeFiles/TT.dir/requires: CMakeFiles/TT.dir/TypeTwo.cpp.o.requires
.PHONY : CMakeFiles/TT.dir/requires

CMakeFiles/TT.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/TT.dir/cmake_clean.cmake
.PHONY : CMakeFiles/TT.dir/clean

CMakeFiles/TT.dir/depend:
	cd /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/cmakedir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0 /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0 /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/limSend_1.0/cmakedir/CMakeFiles/TT.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TT.dir/depend
