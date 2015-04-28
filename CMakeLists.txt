CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.3 HINTS ${DEAL_II_DIR} )

IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

#cmake macros defined in : build_3/share/deal.II/macros/
DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(ALE)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

#ADD_DEFINITIONS(-g)

INCLUDE_DIRECTORIES(include)

#LINK_DIRECTORIES(${ALE_SOURCE_DIR}/lib)

ADD_SUBDIRECTORY(src)
