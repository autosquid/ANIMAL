include_directories(${CMAKE_CURRENT_SOURCE_DIR}/)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/)

include (GenerateExportHeader)

file(GLOB_RECURSE h_s "maths/*.h")
file(GLOB_RECURSE x_s *.cpp)

message (allhs "${h_s}")

add_library(CAPGMath ${LIB_TYPE} ${x_s})

target_link_libraries(CAPGMath)

# GENERATE_EXPORT_HEADER(CAPGMath)

install(TARGETS CAPGMath
  EXPORT CAPGMathTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
)

install(FILES ${h_s} DESTINATION "${INSTALL_INCLUDE_DIR}/CAPGMath/maths")
install(FILES "maths.h" DESTINATION "${INSTALL_INCLUDE_DIR}/CAPGMath")
