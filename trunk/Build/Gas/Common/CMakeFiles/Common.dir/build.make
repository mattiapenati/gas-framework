# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.6

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/penaz/Documents/Polimi/Progetti/GasFramework

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/penaz/Documents/Polimi/Progetti/GasFramework/Build

# Include any dependencies generated for this target.
include Gas/Common/CMakeFiles/Common.dir/depend.make

# Include the progress variables for this target.
include Gas/Common/CMakeFiles/Common.dir/progress.make

# Include the compile flags for this target's objects.
include Gas/Common/CMakeFiles/Common.dir/flags.make

Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o: Gas/Common/CMakeFiles/Common.dir/flags.make
Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o: ../Gas/Common/Exception.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Common.dir/Exception.cpp.o -c /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Exception.cpp

Gas/Common/CMakeFiles/Common.dir/Exception.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Common.dir/Exception.cpp.i"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Exception.cpp > CMakeFiles/Common.dir/Exception.cpp.i

Gas/Common/CMakeFiles/Common.dir/Exception.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Common.dir/Exception.cpp.s"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Exception.cpp -o CMakeFiles/Common.dir/Exception.cpp.s

Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.requires:
.PHONY : Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.requires

Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.provides: Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.requires
	$(MAKE) -f Gas/Common/CMakeFiles/Common.dir/build.make Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.provides.build
.PHONY : Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.provides

Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.provides.build: Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o
.PHONY : Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.provides.build

Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o: Gas/Common/CMakeFiles/Common.dir/flags.make
Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o: ../Gas/Common/Iterable.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Common.dir/Iterable.cpp.o -c /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Iterable.cpp

Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Common.dir/Iterable.cpp.i"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Iterable.cpp > CMakeFiles/Common.dir/Iterable.cpp.i

Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Common.dir/Iterable.cpp.s"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Iterable.cpp -o CMakeFiles/Common.dir/Iterable.cpp.s

Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.requires:
.PHONY : Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.requires

Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.provides: Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.requires
	$(MAKE) -f Gas/Common/CMakeFiles/Common.dir/build.make Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.provides.build
.PHONY : Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.provides

Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.provides.build: Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o
.PHONY : Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.provides.build

Gas/Common/CMakeFiles/Common.dir/Program.cpp.o: Gas/Common/CMakeFiles/Common.dir/flags.make
Gas/Common/CMakeFiles/Common.dir/Program.cpp.o: ../Gas/Common/Program.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Gas/Common/CMakeFiles/Common.dir/Program.cpp.o"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/Common.dir/Program.cpp.o -c /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Program.cpp

Gas/Common/CMakeFiles/Common.dir/Program.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Common.dir/Program.cpp.i"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Program.cpp > CMakeFiles/Common.dir/Program.cpp.i

Gas/Common/CMakeFiles/Common.dir/Program.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Common.dir/Program.cpp.s"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common/Program.cpp -o CMakeFiles/Common.dir/Program.cpp.s

Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.requires:
.PHONY : Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.requires

Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.provides: Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.requires
	$(MAKE) -f Gas/Common/CMakeFiles/Common.dir/build.make Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.provides.build
.PHONY : Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.provides

Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.provides.build: Gas/Common/CMakeFiles/Common.dir/Program.cpp.o
.PHONY : Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.provides.build

# Object files for target Common
Common_OBJECTS = \
"CMakeFiles/Common.dir/Exception.cpp.o" \
"CMakeFiles/Common.dir/Iterable.cpp.o" \
"CMakeFiles/Common.dir/Program.cpp.o"

# External object files for target Common
Common_EXTERNAL_OBJECTS =

Gas/Common/libCommon.a: Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o
Gas/Common/libCommon.a: Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o
Gas/Common/libCommon.a: Gas/Common/CMakeFiles/Common.dir/Program.cpp.o
Gas/Common/libCommon.a: Gas/Common/CMakeFiles/Common.dir/build.make
Gas/Common/libCommon.a: Gas/Common/CMakeFiles/Common.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library libCommon.a"
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && $(CMAKE_COMMAND) -P CMakeFiles/Common.dir/cmake_clean_target.cmake
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Common.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Gas/Common/CMakeFiles/Common.dir/build: Gas/Common/libCommon.a
.PHONY : Gas/Common/CMakeFiles/Common.dir/build

Gas/Common/CMakeFiles/Common.dir/requires: Gas/Common/CMakeFiles/Common.dir/Exception.cpp.o.requires
Gas/Common/CMakeFiles/Common.dir/requires: Gas/Common/CMakeFiles/Common.dir/Iterable.cpp.o.requires
Gas/Common/CMakeFiles/Common.dir/requires: Gas/Common/CMakeFiles/Common.dir/Program.cpp.o.requires
.PHONY : Gas/Common/CMakeFiles/Common.dir/requires

Gas/Common/CMakeFiles/Common.dir/clean:
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common && $(CMAKE_COMMAND) -P CMakeFiles/Common.dir/cmake_clean.cmake
.PHONY : Gas/Common/CMakeFiles/Common.dir/clean

Gas/Common/CMakeFiles/Common.dir/depend:
	cd /home/penaz/Documents/Polimi/Progetti/GasFramework/Build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/penaz/Documents/Polimi/Progetti/GasFramework /home/penaz/Documents/Polimi/Progetti/GasFramework/Gas/Common /home/penaz/Documents/Polimi/Progetti/GasFramework/Build /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common /home/penaz/Documents/Polimi/Progetti/GasFramework/Build/Gas/Common/CMakeFiles/Common.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Gas/Common/CMakeFiles/Common.dir/depend

