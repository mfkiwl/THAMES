cmake_minimum_required(VERSION 3.30)

project (THAMES)

set (THAMES_VERSION 3.0)
set (THAMES_VERSION_MAJOR 3)
set (THAMES_VERSION_MINOR 0)

set (thames_XSDFILES_DIR "${CMAKE_SOURCE_DIR}/src/thameslib/xsdfiles")

# configure a header file to pass some CMake settings to source code
configure_file (
    "${CMAKE_SOURCE_DIR}/src/version.h.in"
    "${CMAKE_SOURCE_DIR}/src/version.h"
    )

configure_file (
    "${CMAKE_SOURCE_DIR}/src/vcctl2thames/vcctl2thames.h.in"
    "${CMAKE_SOURCE_DIR}/src/vcctl2thames/vcctl2thames.h"
    )

option(BUILD_DOC "Build API documentation" ON)
# look for doxygen capability
find_package(Doxygen
              REQUIRED dot
              OPTIONAL_COMPONENTS mscgen dia)

if (DOXYGEN_FOUND)
  # set input and output files
  # set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/src/Doxyfile.in)
  # set(DOXYGEN_OUT ${CMAKE_SOURCE_DIR}/src/Doxyfile)
  # request to configure the file
  #configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT}$ @ONLY)
  doxygen_add_docs(doc
    filesOrDirs ${CMAKE_SOURCE_DIR}/src
    COMMENT "Generating API documentation with Doxygen" VERBATIM
  )

  # message("Doxygen build started")

  # add_custom_target (doc
  #   COMMAND ${Doxygen::doxygen} ${DOXYGEN_OUT}
  #   COMMENT "Generating API documentation with Doxygen" VERBATIM
  # )
else (DOXYGEN_FOUND)
  message("Doxygen must be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

find_program (
  PDFLATEX
  NAMES pdflatex
  HINTS "/Library/TeX/texbin /usr/bin /usr/local/bin"
)

if (PDFLATEX)
  message("pdflatex found... yes")
else()
  message("pdflatex found... no")
endif()

find_library (
    GEMS3K_LIB
    NAMES GEMS3K-static libGEMS3K-static libGEMS3K-static.a GEMS3K-static.a # what to look for
    HINTS "${CMAKE_SOURCE_DIR}/src/Resources/lib" # where to look
)

# check if we found the GEMS3K library
message(STATUS "GEMS3K_LIB: [${GEMS3K_LIB}]")

if (GEMS3K_LIB)
    message ("Found GEMSK3K_LIB: [${GEMS3K_LIB}]")
    set (EXTRA_LIBS ${EXTRA_LIBS} ${GEMS3K_LIB})
    message ("EXTRA_LIBS: " ${EXTRA_LIBS})
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

set (THAMES_INCLUDE_PATH "-I${CMAKE_SOURCE_DIR}/src -I/usr/local/include")

set (CMAKE_CXX_FLAGS "-O2 -std=c++17 -DIPMGEMPLUGIN ${THAMES_INCLUDE_PATH}")

message ("Top level CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message ("Linked libraries are: " ${EXTRA_LIBS})

add_subdirectory (${CMAKE_SOURCE_DIR}/src/thameslib)
add_subdirectory (${CMAKE_SOURCE_DIR}/src/vcctl2thames)
add_subdirectory (${CMAKE_SOURCE_DIR}/src/viz)

file (GLOB SOURCES "${CMAKE_SOURCE_DIR}/src/*.cc")

#set (EXTRA_LIBS ${EXTRA_LIBS} z)
#set (EXTRA_LIBS ${EXTRA_LIBS} m)

add_executable (thames ${SOURCES})

target_link_libraries (thames thameslib ${EXTRA_LIBS})

install (TARGETS thames viz vcctl2thames DESTINATION ${CMAKE_SOURCE_DIR}/bin)
install (DIRECTORY ${CMAKE_SOURCE_DIR}/build/html DESTINATION ${CMAKE_SOURCE_DIR}/doc)
