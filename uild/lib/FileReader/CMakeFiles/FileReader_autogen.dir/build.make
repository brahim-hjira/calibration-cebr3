# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/hjira/calibration/TEMPORAL-main

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/hjira/calibration/TEMPORAL-main/uild

# Utility rule file for FileReader_autogen.

# Include the progress variables for this target.
include lib/FileReader/CMakeFiles/FileReader_autogen.dir/progress.make

lib/FileReader/CMakeFiles/FileReader_autogen:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/hjira/calibration/TEMPORAL-main/uild/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Automatic MOC and UIC for target FileReader"
	cd /home/hjira/calibration/TEMPORAL-main/uild/lib/FileReader && /usr/bin/cmake -E cmake_autogen /home/hjira/calibration/TEMPORAL-main/uild/lib/FileReader/CMakeFiles/FileReader_autogen.dir/AutogenInfo.json ""

FileReader_autogen: lib/FileReader/CMakeFiles/FileReader_autogen
FileReader_autogen: lib/FileReader/CMakeFiles/FileReader_autogen.dir/build.make

.PHONY : FileReader_autogen

# Rule to build all files generated by this target.
lib/FileReader/CMakeFiles/FileReader_autogen.dir/build: FileReader_autogen

.PHONY : lib/FileReader/CMakeFiles/FileReader_autogen.dir/build

lib/FileReader/CMakeFiles/FileReader_autogen.dir/clean:
	cd /home/hjira/calibration/TEMPORAL-main/uild/lib/FileReader && $(CMAKE_COMMAND) -P CMakeFiles/FileReader_autogen.dir/cmake_clean.cmake
.PHONY : lib/FileReader/CMakeFiles/FileReader_autogen.dir/clean

lib/FileReader/CMakeFiles/FileReader_autogen.dir/depend:
	cd /home/hjira/calibration/TEMPORAL-main/uild && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/hjira/calibration/TEMPORAL-main /home/hjira/calibration/TEMPORAL-main/lib/FileReader /home/hjira/calibration/TEMPORAL-main/uild /home/hjira/calibration/TEMPORAL-main/uild/lib/FileReader /home/hjira/calibration/TEMPORAL-main/uild/lib/FileReader/CMakeFiles/FileReader_autogen.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/FileReader/CMakeFiles/FileReader_autogen.dir/depend

