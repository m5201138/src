# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.8

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/m5201138/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/m5201138/src

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/opt/local/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/opt/local/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/m5201138/src/CMakeFiles /Users/m5201138/src/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/m5201138/src/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named pick

# Build rule for target.
pick: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 pick
.PHONY : pick

# fast build rule for target.
pick/fast:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/build
.PHONY : pick/fast

CIsoSurface.o: CIsoSurface.cpp.o

.PHONY : CIsoSurface.o

# target to build an object file
CIsoSurface.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/CIsoSurface.cpp.o
.PHONY : CIsoSurface.cpp.o

CIsoSurface.i: CIsoSurface.cpp.i

.PHONY : CIsoSurface.i

# target to preprocess a source file
CIsoSurface.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/CIsoSurface.cpp.i
.PHONY : CIsoSurface.cpp.i

CIsoSurface.s: CIsoSurface.cpp.s

.PHONY : CIsoSurface.s

# target to generate assembly for a file
CIsoSurface.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/CIsoSurface.cpp.s
.PHONY : CIsoSurface.cpp.s

Camera.o: Camera.cpp.o

.PHONY : Camera.o

# target to build an object file
Camera.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Camera.cpp.o
.PHONY : Camera.cpp.o

Camera.i: Camera.cpp.i

.PHONY : Camera.i

# target to preprocess a source file
Camera.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Camera.cpp.i
.PHONY : Camera.cpp.i

Camera.s: Camera.cpp.s

.PHONY : Camera.s

# target to generate assembly for a file
Camera.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Camera.cpp.s
.PHONY : Camera.cpp.s

Image.o: Image.cpp.o

.PHONY : Image.o

# target to build an object file
Image.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Image.cpp.o
.PHONY : Image.cpp.o

Image.i: Image.cpp.i

.PHONY : Image.i

# target to preprocess a source file
Image.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Image.cpp.i
.PHONY : Image.cpp.i

Image.s: Image.cpp.s

.PHONY : Image.s

# target to generate assembly for a file
Image.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Image.cpp.s
.PHONY : Image.cpp.s

Shader.o: Shader.cpp.o

.PHONY : Shader.o

# target to build an object file
Shader.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Shader.cpp.o
.PHONY : Shader.cpp.o

Shader.i: Shader.cpp.i

.PHONY : Shader.i

# target to preprocess a source file
Shader.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Shader.cpp.i
.PHONY : Shader.cpp.i

Shader.s: Shader.cpp.s

.PHONY : Shader.s

# target to generate assembly for a file
Shader.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Shader.cpp.s
.PHONY : Shader.cpp.s

TriMesh.o: TriMesh.cpp.o

.PHONY : TriMesh.o

# target to build an object file
TriMesh.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/TriMesh.cpp.o
.PHONY : TriMesh.cpp.o

TriMesh.i: TriMesh.cpp.i

.PHONY : TriMesh.i

# target to preprocess a source file
TriMesh.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/TriMesh.cpp.i
.PHONY : TriMesh.cpp.i

TriMesh.s: TriMesh.cpp.s

.PHONY : TriMesh.s

# target to generate assembly for a file
TriMesh.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/TriMesh.cpp.s
.PHONY : TriMesh.cpp.s

Vectors.o: Vectors.cpp.o

.PHONY : Vectors.o

# target to build an object file
Vectors.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Vectors.cpp.o
.PHONY : Vectors.cpp.o

Vectors.i: Vectors.cpp.i

.PHONY : Vectors.i

# target to preprocess a source file
Vectors.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Vectors.cpp.i
.PHONY : Vectors.cpp.i

Vectors.s: Vectors.cpp.s

.PHONY : Vectors.s

# target to generate assembly for a file
Vectors.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Vectors.cpp.s
.PHONY : Vectors.cpp.s

Viewer.o: Viewer.cpp.o

.PHONY : Viewer.o

# target to build an object file
Viewer.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Viewer.cpp.o
.PHONY : Viewer.cpp.o

Viewer.i: Viewer.cpp.i

.PHONY : Viewer.i

# target to preprocess a source file
Viewer.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Viewer.cpp.i
.PHONY : Viewer.cpp.i

Viewer.s: Viewer.cpp.s

.PHONY : Viewer.s

# target to generate assembly for a file
Viewer.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/Viewer.cpp.s
.PHONY : Viewer.cpp.s

main.o: main.cpp.o

.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i

.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s

.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/main.cpp.s
.PHONY : main.cpp.s

types.o: types.cpp.o

.PHONY : types.o

# target to build an object file
types.cpp.o:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/types.cpp.o
.PHONY : types.cpp.o

types.i: types.cpp.i

.PHONY : types.i

# target to preprocess a source file
types.cpp.i:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/types.cpp.i
.PHONY : types.cpp.i

types.s: types.cpp.s

.PHONY : types.s

# target to generate assembly for a file
types.cpp.s:
	$(MAKE) -f CMakeFiles/pick.dir/build.make CMakeFiles/pick.dir/types.cpp.s
.PHONY : types.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... pick"
	@echo "... CIsoSurface.o"
	@echo "... CIsoSurface.i"
	@echo "... CIsoSurface.s"
	@echo "... Camera.o"
	@echo "... Camera.i"
	@echo "... Camera.s"
	@echo "... Image.o"
	@echo "... Image.i"
	@echo "... Image.s"
	@echo "... Shader.o"
	@echo "... Shader.i"
	@echo "... Shader.s"
	@echo "... TriMesh.o"
	@echo "... TriMesh.i"
	@echo "... TriMesh.s"
	@echo "... Vectors.o"
	@echo "... Vectors.i"
	@echo "... Vectors.s"
	@echo "... Viewer.o"
	@echo "... Viewer.i"
	@echo "... Viewer.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... types.o"
	@echo "... types.i"
	@echo "... types.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

