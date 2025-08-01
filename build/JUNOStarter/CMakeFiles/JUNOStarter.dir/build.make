# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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
CMAKE_COMMAND = /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J24.2.0/ExternalLibs/Cmake/3.30.5/bin/cmake

# The command to remove a file.
RM = /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J24.2.0/ExternalLibs/Cmake/3.30.5/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build

# Include any dependencies generated for this target.
include JUNOStarter/CMakeFiles/JUNOStarter.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include JUNOStarter/CMakeFiles/JUNOStarter.dir/compiler_depend.make

# Include the progress variables for this target.
include JUNOStarter/CMakeFiles/JUNOStarter.dir/progress.make

# Include the compile flags for this target's objects.
include JUNOStarter/CMakeFiles/JUNOStarter.dir/flags.make

JUNOStarter/CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o: JUNOStarter/CMakeFiles/JUNOStarter.dir/flags.make
JUNOStarter/CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o: /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/JUNOStarter/src/JUNOStarter.cc
JUNOStarter/CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o: JUNOStarter/CMakeFiles/JUNOStarter.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object JUNOStarter/CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o"
	cd /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/JUNOStarter && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT JUNOStarter/CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o -MF CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o.d -o CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o -c /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/JUNOStarter/src/JUNOStarter.cc

JUNOStarter/CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.i"
	cd /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/JUNOStarter && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/JUNOStarter/src/JUNOStarter.cc > CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.i

JUNOStarter/CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.s"
	cd /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/JUNOStarter && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/JUNOStarter/src/JUNOStarter.cc -o CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.s

# Object files for target JUNOStarter
JUNOStarter_OBJECTS = \
"CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o"

# External object files for target JUNOStarter
JUNOStarter_EXTERNAL_OBJECTS =

lib/libJUNOStarter.so: JUNOStarter/CMakeFiles/JUNOStarter.dir/src/JUNOStarter.cc.o
lib/libJUNOStarter.so: JUNOStarter/CMakeFiles/JUNOStarter.dir/build.make
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libEvtNavigator.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libBufferMemMgrLib.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/sniper/InstallArea/lib64/libRootWriterLib.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libGeometryLib.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libCalibEvent.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libElecEvent.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libRecEvent.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/sniper/InstallArea/lib64/libSniperKernel.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libCore.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libImt.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libRIO.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libNet.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libHist.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libGraf.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libGraf3d.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libGpad.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libROOTDataFrame.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libTree.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libTreePlayer.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libRint.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libPostscript.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libMatrix.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libPhysics.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libMathCore.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libThread.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libMultiProc.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libROOTVecOps.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libCore.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libImt.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libRIO.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libNet.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libHist.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libGraf.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libGraf3d.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libGpad.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libROOTDataFrame.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libTree.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libTreePlayer.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libRint.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libPostscript.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libMatrix.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libPhysics.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libMathCore.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libThread.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libMultiProc.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.4/ExternalLibs/ROOT/6.30.08/lib/libROOTVecOps.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libIdentifier.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J24.2.0/ExternalLibs/Boost/1.85.0/lib/libboost_system.so.1.85.0
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J24.2.0/ExternalLibs/Boost/1.85.0/lib/libboost_python311.so.1.85.0
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/Python/3.11.10/lib/libpython3.11.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libDataPathHelper.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libGeom.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libHist.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libPhysics.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libMatrix.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libGenVector.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libMathCore.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libImt.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libMultiProc.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libNet.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libRIO.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libThread.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libContext.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libCLHEPDict.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libEDMUtil.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/junosw/InstallArea/lib64/libBaseEvent.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/ROOT/6.30.08/lib/libCore.so
lib/libJUNOStarter.so: /cvmfs/juno.ihep.ac.cn/el9_amd64_gcc11/Release/J25.1.6/ExternalLibs/CLHEP/2.4.7.1/lib/libCLHEP.so
lib/libJUNOStarter.so: JUNOStarter/CMakeFiles/JUNOStarter.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared module ../lib/libJUNOStarter.so"
	cd /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/JUNOStarter && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/JUNOStarter.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
JUNOStarter/CMakeFiles/JUNOStarter.dir/build: lib/libJUNOStarter.so
.PHONY : JUNOStarter/CMakeFiles/JUNOStarter.dir/build

JUNOStarter/CMakeFiles/JUNOStarter.dir/clean:
	cd /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/JUNOStarter && $(CMAKE_COMMAND) -P CMakeFiles/JUNOStarter.dir/cmake_clean.cmake
.PHONY : JUNOStarter/CMakeFiles/JUNOStarter.dir/clean

JUNOStarter/CMakeFiles/JUNOStarter.dir/depend:
	cd /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/JUNOStarter /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/JUNOStarter /datafs/users/wujxy/waterphase_analysis/realdata_analysis/juno_analysis/juno-starter/build/JUNOStarter/CMakeFiles/JUNOStarter.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : JUNOStarter/CMakeFiles/JUNOStarter.dir/depend

