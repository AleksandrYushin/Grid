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
CMAKE_SOURCE_DIR = /home/alexandro/Programms/VS_code_progect/Setca

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/alexandro/Programms/VS_code_progect/Setca/build

# Include any dependencies generated for this target.
include CMakeFiles/cub.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/cub.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/cub.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cub.dir/flags.make

CMakeFiles/cub.dir/cub.cpp.o: CMakeFiles/cub.dir/flags.make
CMakeFiles/cub.dir/cub.cpp.o: ../cub.cpp
CMakeFiles/cub.dir/cub.cpp.o: CMakeFiles/cub.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/alexandro/Programms/VS_code_progect/Setca/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cub.dir/cub.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/cub.dir/cub.cpp.o -MF CMakeFiles/cub.dir/cub.cpp.o.d -o CMakeFiles/cub.dir/cub.cpp.o -c /home/alexandro/Programms/VS_code_progect/Setca/cub.cpp

CMakeFiles/cub.dir/cub.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cub.dir/cub.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/alexandro/Programms/VS_code_progect/Setca/cub.cpp > CMakeFiles/cub.dir/cub.cpp.i

CMakeFiles/cub.dir/cub.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cub.dir/cub.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/alexandro/Programms/VS_code_progect/Setca/cub.cpp -o CMakeFiles/cub.dir/cub.cpp.s

# Object files for target cub
cub_OBJECTS = \
"CMakeFiles/cub.dir/cub.cpp.o"

# External object files for target cub
cub_EXTERNAL_OBJECTS =

cub: CMakeFiles/cub.dir/cub.cpp.o
cub: CMakeFiles/cub.dir/build.make
cub: /usr/lib/x86_64-linux-gnu/libgmsh.so
cub: CMakeFiles/cub.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/alexandro/Programms/VS_code_progect/Setca/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable cub"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cub.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cub.dir/build: cub
.PHONY : CMakeFiles/cub.dir/build

CMakeFiles/cub.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cub.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cub.dir/clean

CMakeFiles/cub.dir/depend:
	cd /home/alexandro/Programms/VS_code_progect/Setca/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/alexandro/Programms/VS_code_progect/Setca /home/alexandro/Programms/VS_code_progect/Setca /home/alexandro/Programms/VS_code_progect/Setca/build /home/alexandro/Programms/VS_code_progect/Setca/build /home/alexandro/Programms/VS_code_progect/Setca/build/CMakeFiles/cub.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cub.dir/depend

