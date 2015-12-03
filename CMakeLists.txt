cmake_minimum_required(VERSION 2.6)
project(CAPGMath)


# Offer the user the choice of overriding the installation directories
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR CMake)
else()
  set(DEF_INSTALL_CMAKE_DIR lib/CMake/CAPGMath)
endif()

set(INSTALL_CMAKE_DIR ${DEF_INSTALL_CMAKE_DIR} CACHE PATH "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INCLUDE CMAKE)
  set(var INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" ON)
set (CMAKE_DEBUG_POSTFIX "d")

IF (BUILD_SHARED_LIBS)
  # User wants to build Dynamic Libraries, so change the LIB_TYPE variable to CMake keyword 'SHARED'
  SET (LIB_TYPE SHARED)
ELSE(BUILD_SHARED_LIBS)
  SET(LIB_TPYE STATIC)
ENDIF (BUILD_SHARED_LIBS)


#add source directory
add_subdirectory(src)

# The interesting stuff goes here
# ===============================

# Add all targets to the build-tree export set
export(TARGETS CAPGMath FILE "${PROJECT_BINARY_DIR}/CAPGMathTargets.cmake")

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE CAPGMath)

file(RELATIVE_PATH REL_INCLUDE_DIR "${INSTALL_CMAKE_DIR}" "${INSTALL_INCLUDE_DIR}")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
configure_file(CAPGMathConfig.cmake.in "${PROJECT_BINARY_DIR}/CAPGMathConfig.cmake" @ONLY)

# ... for the install tree
set(CONF_INCLUDE_DIRS "\${CAPGMath_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(CAPGMathConfig.cmake.in "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CAPGMathConfig.cmake" @ONLY)

install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CAPGMathConfig.cmake" DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)

# Install the export set for use with the install-tree
install(EXPORT CAPGMathTargets DESTINATION "${INSTALL_CMAKE_DIR}" COMPONENT dev)
