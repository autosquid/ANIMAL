include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/)
include (GenerateExportHeader)

file(GLOB_RECURSE h_s *.h)
file(GLOB_RECURSE x_s *.cpp)


OPTION (BUILD_SHARED_LIBS "Build Shared Libraries" ON)
set (CMAKE_DEBUG_POSTFIX "d")

IF (BUILD_SHARED_LIBS)
  # User wants to build Dynamic Libraries, so change the LIB_TYPE variable to CMake keyword 'SHARED'
  SET (LIB_TYPE SHARED)
ELSE(BUILD_SHARED_LIBS)
	SET(LIB_TPYE STATIC)
ENDIF (BUILD_SHARED_LIBS)

add_library(CAPGMath ${LIB_TYPE} ${h_s} ${x_s}) # ${moc_s} ${u_h_s})

GENERATE_EXPORT_HEADER(CAPGMath)
target_link_libraries(CAPGMath) # ${QT_LIBRARIES} ${OpenCV_LIBS} ${GLEW_LIBRARY})

