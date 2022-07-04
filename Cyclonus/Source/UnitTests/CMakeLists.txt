ENABLE_TESTING()

ADD_EXECUTABLE(Check
  ProtoFiles/PhysicalConstants/Fundamental_tests.cpp
)

TARGET_LINK_LIBRARIES(Check
  ProtoBinding
  gtest_main
  )

INCLUDE(GoogleTest)
gtest_discover_tests(Check)

SET_TARGET_PROPERTIES(Check
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/")