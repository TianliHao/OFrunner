# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/tlh/Dev/OpenMesh-8.0

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/tlh/Dev/OpenMesh-8.0/BUILD

# Utility rule file for DecimaterGui_autogen.

# Include the progress variables for this target.
include src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/progress.make

src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/tlh/Dev/OpenMesh-8.0/BUILD/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Automatic MOC for target DecimaterGui"
	cd /home/tlh/Dev/OpenMesh-8.0/BUILD/src/OpenMesh/Apps/Decimating/DecimaterGui && /usr/bin/cmake -E cmake_autogen /home/tlh/Dev/OpenMesh-8.0/BUILD/src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir Release

DecimaterGui_autogen: src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen
DecimaterGui_autogen: src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/build.make

.PHONY : DecimaterGui_autogen

# Rule to build all files generated by this target.
src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/build: DecimaterGui_autogen

.PHONY : src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/build

src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/clean:
	cd /home/tlh/Dev/OpenMesh-8.0/BUILD/src/OpenMesh/Apps/Decimating/DecimaterGui && $(CMAKE_COMMAND) -P CMakeFiles/DecimaterGui_autogen.dir/cmake_clean.cmake
.PHONY : src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/clean

src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/depend:
	cd /home/tlh/Dev/OpenMesh-8.0/BUILD && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/tlh/Dev/OpenMesh-8.0 /home/tlh/Dev/OpenMesh-8.0/src/OpenMesh/Apps/Decimating/DecimaterGui /home/tlh/Dev/OpenMesh-8.0/BUILD /home/tlh/Dev/OpenMesh-8.0/BUILD/src/OpenMesh/Apps/Decimating/DecimaterGui /home/tlh/Dev/OpenMesh-8.0/BUILD/src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/OpenMesh/Apps/Decimating/DecimaterGui/CMakeFiles/DecimaterGui_autogen.dir/depend

