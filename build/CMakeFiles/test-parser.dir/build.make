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
CMAKE_SOURCE_DIR = /home/othman/Equipe_0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/othman/Equipe_0/build

# Include any dependencies generated for this target.
include CMakeFiles/test-parser.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/test-parser.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/test-parser.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/test-parser.dir/flags.make

CMakeFiles/test-parser.dir/src/test_parser.cpp.o: CMakeFiles/test-parser.dir/flags.make
CMakeFiles/test-parser.dir/src/test_parser.cpp.o: ../src/test_parser.cpp
CMakeFiles/test-parser.dir/src/test_parser.cpp.o: CMakeFiles/test-parser.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/othman/Equipe_0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/test-parser.dir/src/test_parser.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test-parser.dir/src/test_parser.cpp.o -MF CMakeFiles/test-parser.dir/src/test_parser.cpp.o.d -o CMakeFiles/test-parser.dir/src/test_parser.cpp.o -c /home/othman/Equipe_0/src/test_parser.cpp

CMakeFiles/test-parser.dir/src/test_parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-parser.dir/src/test_parser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/othman/Equipe_0/src/test_parser.cpp > CMakeFiles/test-parser.dir/src/test_parser.cpp.i

CMakeFiles/test-parser.dir/src/test_parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-parser.dir/src/test_parser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/othman/Equipe_0/src/test_parser.cpp -o CMakeFiles/test-parser.dir/src/test_parser.cpp.s

CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o: CMakeFiles/test-parser.dir/flags.make
CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o: ../src/3rdparty/jlparser/parser.cpp
CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o: CMakeFiles/test-parser.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/othman/Equipe_0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o -MF CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o.d -o CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o -c /home/othman/Equipe_0/src/3rdparty/jlparser/parser.cpp

CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/othman/Equipe_0/src/3rdparty/jlparser/parser.cpp > CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.i

CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/othman/Equipe_0/src/3rdparty/jlparser/parser.cpp -o CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.s

# Object files for target test-parser
test__parser_OBJECTS = \
"CMakeFiles/test-parser.dir/src/test_parser.cpp.o" \
"CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o"

# External object files for target test-parser
test__parser_EXTERNAL_OBJECTS =

test-parser: CMakeFiles/test-parser.dir/src/test_parser.cpp.o
test-parser: CMakeFiles/test-parser.dir/src/3rdparty/jlparser/parser.cpp.o
test-parser: CMakeFiles/test-parser.dir/build.make
test-parser: /home/othman/pnl-1.11.0/build/lib/libpnl.so
test-parser: /usr/lib/x86_64-linux-gnu/libopenblas.so
test-parser: /usr/lib/x86_64-linux-gnu/libopenblas.so
test-parser: /home/othman/anaconda3/lib/libmpi.so
test-parser: CMakeFiles/test-parser.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/othman/Equipe_0/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable test-parser"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/test-parser.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/test-parser.dir/build: test-parser
.PHONY : CMakeFiles/test-parser.dir/build

CMakeFiles/test-parser.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/test-parser.dir/cmake_clean.cmake
.PHONY : CMakeFiles/test-parser.dir/clean

CMakeFiles/test-parser.dir/depend:
	cd /home/othman/Equipe_0/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/othman/Equipe_0 /home/othman/Equipe_0 /home/othman/Equipe_0/build /home/othman/Equipe_0/build /home/othman/Equipe_0/build/CMakeFiles/test-parser.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/test-parser.dir/depend

