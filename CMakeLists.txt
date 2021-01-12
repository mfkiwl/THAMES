cmake_minimum_required(VERSION 3.0)

project (thames)

set (THAMES_VERSION 2.5)
set (THAMES_VERSION_MAJOR 2)
set (THAMES_VERSION_MINOR 5.1)

set (thames_XSDFILES_DIR "${CMAKE_SOURCE_DIR}/src/thameslib/xsdfiles")

# configure a header file to pass some CMake settings to source code
configure_file (
    "${CMAKE_SOURCE_DIR}/src/version.h.in"
    "${CMAKE_SOURCE_DIR}/src/version.h"
    )

configure_file (
    "${CMAKE_SOURCE_DIR}/src/thameslib/myconfig.h.in"
    "${CMAKE_SOURCE_DIR}/src/thameslib/myconfig.h"
    )

configure_file (
    "${CMAKE_SOURCE_DIR}/src/vcctl2thames/vcctl2thames.h.in"
    "${CMAKE_SOURCE_DIR}/src/vcctl2thames/vcctl2thames.h"
    )

# look for doxygen capability
find_package (Doxygen)
if (DOXYGEN_FOUND)
    configure_file ("${CMAKE_SOURCE_DIR}/Doxyfile.in" "${CMAKE_BINARY_DIR}/Doxyfile" @ONLY)
    add_custom_target (doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif (DOXYGEN_FOUND)

find_package (LibXml2)
if (LIBXML2_FOUND)
    include_directories(${LIBXML2_INCLUDE_DIR})
    message ("Found libxml2 include files at " ${LIBXML2_INCLUDE_DIR})
    message ("libxml2 libraries are " ${LIBXML2_LIBRARIES})
    set (EXTRA_LIBS ${EXTRA_LIBS} ${LIBXML2_LIBRARIES})
    set (XML_INCLUDE_PATH "-I${LIBXML2_INCLUDE_DIR}")
else (LIBXML2_FOUND)
    message (FATAL_ERROR "Library libxml not installed on system.  Please install it first.")
endif (LIBXML2_FOUND)

find_library (
    GEMS3K_LIB
    NAMES GEMS3K-static libGEMS3K-static libGEMS3K-static.a GEMS3K-static.a # what to look for
    HINTS "/usr/local/lib" # where to look
)

# check if we found the GEMS3K library
message(STATUS "GEMS3K_LIB: [${GEMS3K_LIB}]")

if (GEMS3K_LIB)
    message ("Found GEMSK3K_LIB: [${GEMS3K_LIB}]")
    set (EXTRA_LIBS ${EXTRA_LIBS} ${GEMS3K_LIB})
else ()
    message (FATAL_ERROR "Did not find lib GEMS3K-static")
endif ()

find_library (
    MATH_LIB
    NAMES m libm # what to look for
    HINTS "/usr/lib" # where to look
)

message(STATUS "MATH_LIB: [${MATH_LIB}]")

if (MATH_LIB)
    message ("Found MATH_LIB: [${MATH_LIB}]")
    set (EXTRA_LIBS ${EXTRA_LIBS} ${MATH_LIB})
else ()
    message (FATAL_ERROR "Did not find MATH lib")
endif ()

find_library (
    Z_LIB
    NAMES z libz # what to look for
    HINTS "/usr/lib" # where to look
)

message(STATUS "Z_LIB: [${Z_LIB}]")

if (Z_LIB)
    message ("Found Z_LIB: [${Z_LIB}]")
    set (EXTRA_LIBS ${EXTRA_LIBS} ${Z_LIB})
else ()
    message (FATAL_ERROR "Did not find Z lib")
endif ()

set (THAMES_INCLUDE_PATH "-I${CMAKE_SOURCE_DIR}/src -I/usr/local/include")

set (CMAKE_CXX_FLAGS "-O2 -std=c++17 -DIPMGEMPLUGIN ${THAMES_INCLUDE_PATH} ${XML_INCLUDE_PATH}")

message ("Top level CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

add_subdirectory (${CMAKE_SOURCE_DIR}/src/thameslib)
add_subdirectory (${CMAKE_SOURCE_DIR}/src/vcctl2thames)

file (GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cc")

#set (EXTRA_LIBS ${EXTRA_LIBS} z)
#set (EXTRA_LIBS ${EXTRA_LIBS} m)
#set (EXTRA_LIBS ${EXTRA_LIBS} xml)

add_executable (thames ${SOURCES})

target_link_libraries (thames thameslib ${EXTRA_LIBS})

install (TARGETS thames DESTINATION ${CMAKE_SOURCE_DIR}/bin)
