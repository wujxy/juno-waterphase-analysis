
cmake_minimum_required(VERSION 3.20)
project(AnalysisHelloWorld)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/InstallArea)
set(CMAKE_BUILD_TYPE Debug)

find_package(junosw REQUIRED CONFIG)
include(PKG)

add_subdirectory(JUNOStarter)
add_subdirectory(JUNOStarterESD)

include(JUNOSetupScripts)
include(JUNOInstallOfflineSetup)