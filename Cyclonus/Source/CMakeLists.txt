CMAKE_MINIMUM_REQUIRED (VERSION 3.14)

# Language of the project is CXX (for C++)
PROJECT(CYCLONUS LANGUAGES CXX)

# BEFORE prepends the directories specified onto current list of directories.
# cf. https://cmake.org/cmake/help/v3.0/command/include_directories.html
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_SOURCE_DIR})

# Include will load and run CMake code of this module, "CheckCXXCompilerFlag."
INCLUDE(CheckCXXCompilerFlag)

# Compiler Flags

# cf. https://cmake.org/cmake/help/latest/command/add_compile_options.html
# Add options to the compilation of source files.
# cf. https://stackoverflow.com/questions/29254877/meaning-of-g-flags-wall-w-werror
# Wall - enable a set of warnings, not all of them.
# extra - enable extra warning.
# pedantic - don't use extensions beyond the standard.
# Werror - every warning treated as an error.
ADD_COMPILE_OPTIONS(-Wall -Wextra -pedantic -Werror)

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

# cf. https://github.com/shaochuan/cmake-protobuf-example/blob/master/proto/CMakeLists.txt
# cf. https://cmake.org/cmake/help/latest/module/FindProtobuf.html
# FindProtobuf - locate and configure Google Protocol Buffers library.
#INCLUDE(FindProtobuf)
#FIND_PACKAGE(Protobuf REQUIRED)
#INCLUDE_DIRECTORIES(${Protobuf_INCLUDE_DIRS})

FIND_PROGRAM(CCACHE_FOUND ccache)

IF(CCACHE_FOUND)
  SET(CMAKE_CXX_COMPILER_LAUNCHER ccache)
ENDIF()

# Calls a CMake module that first, searches the file system for the appropriate
# threads package for this platform, and then sets CMAKE_THREAD_LIBS_INIT
# variable (and other variables as well).
FIND_PACKAGE(Threads REQUIRED)

# Declare dependency on GoogleTest which is downloaded from GitHub.
# cf. https://google.github.io/googletest/quickstart-cmake.html
INCLUDE(FetchContent)
FetchContent_Declare(googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

FetchContent_MakeAvailable(googletest)

# Add_subdirectory call makes use of libraries.
# Add_subdirectory adds the source CMakeLists.txt and code files in the
# directory.
#ADD_SUBDIRECTORY(ProtoFiles)
ADD_SUBDIRECTORY(UnitTests)
