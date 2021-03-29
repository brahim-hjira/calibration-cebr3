######################################################################
# @author      : lopezmar (lopezmar@$HOSTNAME)
# @modified by : o.allegrini (o.allegrini@ipnl.in2p3.fr)
# @file        : CMakeLists
# @created     : jeudi oct. 15, 2020 16:45:17 CEST
######################################################################

cmake_minimum_required(VERSION 3.5)
project(Test)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


#####################################################################
#
#                           ROOT
#
#####################################################################

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)

find_package(ROOT REQUIRED COMPONENTS Spectrum)
include(${ROOT_USE_FILE})


#####################################################################
#
#                           DAMAVAN
#
#####################################################################

add_subdirectory(lib/FileReader)
include_directories(/home/hjira/calibration/TEMPORAL-main/lib/FileReader
/home/hjira/calibration/TEMPORAL-main/lib/FileReader/lib/CommonQuadratorHeader/510
/home/hjira/calibration/TEMPORAL-main/functions
${ROOT_INCLUDE_DIRS})


add_definitions(-DNO_YAML ${ROOT_CXX_FLAGS})

add_executable(Test main.cpp)
add_executable(TestR mainR.cpp)

target_link_libraries(Test FileReader ${ROOT_LIBRARIES})
target_link_libraries(TestR FileReader ${ROOT_LIBRARIES})
