cmake_minimum_required(VERSION 2.6)
project(CAPGMath)



# add Litit Modules
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/misc/cmake)
#include(CAPGCMakeConfig)

# set output place for executable and lib and static lib

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(ARCHIVE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/staticlib)

#add source directory 
add_subdirectory(src)

