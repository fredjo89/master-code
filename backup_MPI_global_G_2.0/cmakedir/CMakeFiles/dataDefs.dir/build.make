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
CMAKE_SOURCE_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/cmakedir

# Include any dependencies generated for this target.
include CMakeFiles/dataDefs.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/dataDefs.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/dataDefs.dir/flags.make

CMakeFiles/dataDefs.dir/data_defs.cpp.o: CMakeFiles/dataDefs.dir/flags.make
CMakeFiles/dataDefs.dir/data_defs.cpp.o: ../data_defs.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/cmakedir/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/dataDefs.dir/data_defs.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/dataDefs.dir/data_defs.cpp.o -c /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/data_defs.cpp

CMakeFiles/dataDefs.dir/data_defs.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/dataDefs.dir/data_defs.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/data_defs.cpp > CMakeFiles/dataDefs.dir/data_defs.cpp.i

CMakeFiles/dataDefs.dir/data_defs.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/dataDefs.dir/data_defs.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/data_defs.cpp -o CMakeFiles/dataDefs.dir/data_defs.cpp.s

CMakeFiles/dataDefs.dir/data_defs.cpp.o.requires:
.PHONY : CMakeFiles/dataDefs.dir/data_defs.cpp.o.requires

CMakeFiles/dataDefs.dir/data_defs.cpp.o.provides: CMakeFiles/dataDefs.dir/data_defs.cpp.o.requires
	$(MAKE) -f CMakeFiles/dataDefs.dir/build.make CMakeFiles/dataDefs.dir/data_defs.cpp.o.provides.build
.PHONY : CMakeFiles/dataDefs.dir/data_defs.cpp.o.provides

CMakeFiles/dataDefs.dir/data_defs.cpp.o.provides.build: CMakeFiles/dataDefs.dir/data_defs.cpp.o

# Object files for target dataDefs
dataDefs_OBJECTS = \
"CMakeFiles/dataDefs.dir/data_defs.cpp.o"

# External object files for target dataDefs
dataDefs_EXTERNAL_OBJECTS =

libdataDefs.a: CMakeFiles/dataDefs.dir/data_defs.cpp.o
libdataDefs.a: CMakeFiles/dataDefs.dir/build.make
libdataDefs.a: CMakeFiles/dataDefs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libdataDefs.a"
	$(CMAKE_COMMAND) -P CMakeFiles/dataDefs.dir/cmake_clean_target.cmake
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/dataDefs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/dataDefs.dir/build: libdataDefs.a
.PHONY : CMakeFiles/dataDefs.dir/build

CMakeFiles/dataDefs.dir/requires: CMakeFiles/dataDefs.dir/data_defs.cpp.o.requires
.PHONY : CMakeFiles/dataDefs.dir/requires

CMakeFiles/dataDefs.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/dataDefs.dir/cmake_clean.cmake
.PHONY : CMakeFiles/dataDefs.dir/clean

CMakeFiles/dataDefs.dir/depend:
	cd /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/cmakedir && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0 /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0 /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/cmakedir /home/shomeb/f/fredjoha/Desktop/master-code/backup_MPI_global_G_2.0/cmakedir/CMakeFiles/dataDefs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/dataDefs.dir/depend

