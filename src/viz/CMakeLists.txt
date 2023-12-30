cmake_minimum_required(VERSION 3.0)

project (Viz_project)

add_executable (viz ${CMAKE_SOURCE_DIR}/src/viz/viz.cc)

#find_package (LibXml2)
#if (LIBXML2_FOUND)
#    include_directories(${LIBXML2_INCLUDE_DIR})
#    message ("Found libxml2 include files at " ${LIBXML2_INCLUDE_DIR})
#    message ("libxml2 libraries are " ${LIBXML2_LIBRARIES})
#    set (EXTRA_LIBS ${EXTRA_LIBS} ${LIBXML2_LIBRARIES})
#    set (XML_INCLUDE_PATH "-I${LIBXML2_INCLUDE_DIR}")
#else (LIBXML2_FOUND)
#    message (FATAL_ERROR "Library libxml not installed on system.  Please install it first.")
#endif (LIBXML2_FOUND)

#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${XML_INCLUDE_PATH}")

set (CMAKE_CXX_FLAGS "-O2 -std=c++17 -I${CMAKE_SOURCE_DIR}/src/viz")

message ("viz directory, CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})

install (TARGETS viz DESTINATION ${CMAKE_SOURCE_DIR}/bin)
