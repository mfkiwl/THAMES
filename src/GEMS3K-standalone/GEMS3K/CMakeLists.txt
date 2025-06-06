# Recursively collect all header files from the current directory
file(GLOB HEADER_FILES *.h
  ../../thames.h
  ../../thameslib/*.h)

# Recursively collect all source files from the current directory
file(GLOB SOURCE_FILES *.cpp
  ../../thames.cc
  ../../thameslib/*.cc)

# Recursively collect all source files from the current directory
#file(GLOB SIMDJSON_FILES ${CMAKE_SOURCE_DIR}/simdjson/simdjson.cpp)

# Compile the C++ source code into object files
add_library(GEMS3K_OBJECT OBJECT ${HEADER_FILES} ${SOURCE_FILES}  )

# Compile the static and shared library
add_library(GEMS3K_STATIC STATIC $<TARGET_OBJECTS:GEMS3K_OBJECT>)
add_library(GEMS3K_SHARED SHARED $<TARGET_OBJECTS:GEMS3K_OBJECT>)

# Set name of the static and shared library
set_target_properties(GEMS3K_STATIC PROPERTIES OUTPUT_NAME GEMS3K-static)
set_target_properties(GEMS3K_SHARED PROPERTIES OUTPUT_NAME GEMS3K)

# Create an install target for the static and shared libraries
install(TARGETS GEMS3K_STATIC DESTINATION "lib" COMPONENT libraries)
install(TARGETS GEMS3K_SHARED DESTINATION "lib" COMPONENT libraries)

# Create an install target for the header files
install(DIRECTORY ${GEMS3K_SOURCE_DIR}
    DESTINATION "include" COMPONENT headers
    FILES_MATCHING PATTERN "*.h")

