# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sagarmalik/Documents/git/FDTD/1d

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sagarmalik/Documents/git/FDTD/1d/build

# Include any dependencies generated for this target.
include CMakeFiles/fdtd_1d_simulation.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/fdtd_1d_simulation.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/fdtd_1d_simulation.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fdtd_1d_simulation.dir/flags.make

CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o: CMakeFiles/fdtd_1d_simulation.dir/flags.make
CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o: ../src/main.cpp
CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o: CMakeFiles/fdtd_1d_simulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sagarmalik/Documents/git/FDTD/1d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o -MF CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o.d -o CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o -c /home/sagarmalik/Documents/git/FDTD/1d/src/main.cpp

CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sagarmalik/Documents/git/FDTD/1d/src/main.cpp > CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.i

CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sagarmalik/Documents/git/FDTD/1d/src/main.cpp -o CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.s

CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o: CMakeFiles/fdtd_1d_simulation.dir/flags.make
CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o: ../src/FDTD_1D.cpp
CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o: CMakeFiles/fdtd_1d_simulation.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sagarmalik/Documents/git/FDTD/1d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o -MF CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o.d -o CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o -c /home/sagarmalik/Documents/git/FDTD/1d/src/FDTD_1D.cpp

CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sagarmalik/Documents/git/FDTD/1d/src/FDTD_1D.cpp > CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.i

CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sagarmalik/Documents/git/FDTD/1d/src/FDTD_1D.cpp -o CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.s

# Object files for target fdtd_1d_simulation
fdtd_1d_simulation_OBJECTS = \
"CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o" \
"CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o"

# External object files for target fdtd_1d_simulation
fdtd_1d_simulation_EXTERNAL_OBJECTS =

fdtd_1d_simulation: CMakeFiles/fdtd_1d_simulation.dir/src/main.cpp.o
fdtd_1d_simulation: CMakeFiles/fdtd_1d_simulation.dir/src/FDTD_1D.cpp.o
fdtd_1d_simulation: CMakeFiles/fdtd_1d_simulation.dir/build.make
fdtd_1d_simulation: CMakeFiles/fdtd_1d_simulation.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sagarmalik/Documents/git/FDTD/1d/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable fdtd_1d_simulation"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fdtd_1d_simulation.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fdtd_1d_simulation.dir/build: fdtd_1d_simulation
.PHONY : CMakeFiles/fdtd_1d_simulation.dir/build

CMakeFiles/fdtd_1d_simulation.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fdtd_1d_simulation.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fdtd_1d_simulation.dir/clean

CMakeFiles/fdtd_1d_simulation.dir/depend:
	cd /home/sagarmalik/Documents/git/FDTD/1d/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sagarmalik/Documents/git/FDTD/1d /home/sagarmalik/Documents/git/FDTD/1d /home/sagarmalik/Documents/git/FDTD/1d/build /home/sagarmalik/Documents/git/FDTD/1d/build /home/sagarmalik/Documents/git/FDTD/1d/build/CMakeFiles/fdtd_1d_simulation.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fdtd_1d_simulation.dir/depend

