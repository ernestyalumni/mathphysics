CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

PROJECT(MANIFOLDS LANGUAGES CXX)

INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR})

# cf. https://cmake.org/cmake/help/v3.12/command/list.html
# List operations, where for Modification,
# LIST(APPEND <list> [<element>...])
# Appends elements to the list.
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/../CMake) 

# Include will load and run CMake code of this module, "CheckCXXCompilerFlag."
INCLUDE(CheckCXXCompilerFlag)

SET(CMAKE_CXX_STANDARD_17)
# Require the required C++ standard from the compiler.
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

FIND_PROGRAM(CCACHE_FOUND ccache)

IF(CCACHE_FOUND)
  SET(CMAKE_CXX_COMPILER_LAUNCHER ccache)
ENDIF()

SET(BOOST_MIN_VERSION "1.61.0")

IF (ALLOW_NEWER_BOOST)
  FIND_PACKAGE(Boost ${BOOST_MIN_VERSION} REQUIRED)
ELSE()
  FIND_PACKAGE(Boost)
ENDIF()

# Include directories so that header files can be found.
INCLUDE_DIRECTORIES(BOOST_INCLUDES ${Boost_INCLUDE_DIRS})

# Add_subdirectory call makes use of libraries.
# Add_subdirectory adds the source CMakeLists.txt and code files in the
# directory.
ADD_SUBDIRECTORY(Algebra)
ADD_SUBDIRECTORY(UnitTests)
