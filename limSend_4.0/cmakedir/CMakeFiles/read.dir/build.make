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
CMAKE_SOURCE_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/cmakedir

# Include any dependencies generated for this target.
include CMakeFiles/read.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/read.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/read.dir/flags.make

CMakeFiles/read.dir/readFromFile.cpp.o: CMakeFiles/read.dir/flags.make
CMakeFiles/read.dir/readFromFile.cpp.o: ../readFromFile.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/cmakedir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/read.dir/readFromFile.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/read.dir/readFromFile.cpp.o -c /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/readFromFile.cpp

CMakeFiles/read.dir/readFromFile.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/read.dir/readFromFile.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/readFromFile.cpp > CMakeFiles/read.dir/readFromFile.cpp.i

CMakeFiles/read.dir/readFromFile.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/read.dir/readFromFile.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/readFromFile.cpp -o CMakeFiles/read.dir/readFromFile.cpp.s

CMakeFiles/read.dir/readFromFile.cpp.o.requires:
.PHONY : CMakeFiles/read.dir/readFromFile.cpp.o.requires

CMakeFiles/read.dir/readFromFile.cpp.o.provides: CMakeFiles/read.dir/readFromFile.cpp.o.requires
	$(MAKE) -f CMakeFiles/read.dir/build.make CMakeFiles/read.dir/readFromFile.cpp.o.provides.build
.PHONY : CMakeFiles/read.dir/readFromFile.cpp.o.provides

CMakeFiles/read.dir/readFromFile.cpp.o.provides.build: CMakeFiles/read.dir/readFromFile.cpp.o

# Object files for target read
read_OBJECTS = \
"CMakeFiles/read.dir/readFromFile.cpp.o"

# External object files for target read
read_EXTERNAL_OBJECTS =

libread.a: CMakeFiles/read.dir/readFromFile.cpp.o
libread.a: CMakeFiles/read.dir/build.make
libread.a: CMakeFiles/read.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libread.a"
	$(CMAKE_COMMAND) -P CMakeFiles/read.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/read.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/read.dir/build: libread.a
.PHONY : CMakeFiles/read.dir/build

CMakeFiles/read.dir/requires: CMakeFiles/read.dir/readFromFile.cpp.o.requires
.PHONY : CMakeFiles/read.dir/requires

CMakeFiles/read.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/read.dir/cmake_clean.cmake
.PHONY : CMakeFiles/read.dir/clean

CMakeFiles/read.dir/depend:
	cd /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/cmakedir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0 /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0 /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/limSend_4.0/cmakedir/CMakeFiles/read.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/read.dir/depend

