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

# cf. https://stackoverflow.com/questions/25451254/getting-cmake-check-cxx-compiler-flag-to-work

CHECK_CXX_COMPILER_FLAG(-std=c++20 COMPILER_SUPPORTS_CXX20)

# Check if Compiler supports the Working draft for ISO C++ 2020 standard.
CHECK_CXX_COMPILER_FLAG(-std=c++2a COMPILER_SUPPORTS_CXX2A)

IF(COMPILER_SUPPORTS_CXX20)
  SET(CMAKE_CXX_STANDARD "20") # C++20
  SET(CMAKE_CXX_STANDARD_REQUIRED ON) # ... is required
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
ELSEIF(COMPILER_SUPPORTS_CXX2A)
  SET(CMAKE_CXX_STANDARD "20") # C++20, but with c++2a
  MESSAGE("Setting to the Working draft of C++20: " ${COMPILER_SUPPORTS_CXX2A})
  SET(CMAKE_CXX_STANDARD_REQUIRED ON) # ... is required
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++2a")
ELSE()
  SET(CMAKE_CXX_STANDARD "17")
  SET(CMAKE_CXX_STANDARD_REQUIRED ON)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
ENDIF()

MESSAGE("CMAKE_CXX_STANDARD: " ${CMAKE_CXX_STANDARD})


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
ADD_SUBDIRECTORY(Utilities)
