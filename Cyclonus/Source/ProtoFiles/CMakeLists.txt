# Expected: mathphysics/Cyclonus/Source/../ProtoFiles
SET(PROTOBUF_INPUT_FILES_BASE_PATH ${CMAKE_SOURCE_DIR}/../ProtoFiles)

SET(PROTOBUF_OUTPUT_FILES_BASE_PATH ${CMAKE_SOURCE_DIR}/ProtoFiles)

SET(PROTOBUF_OUTPUT_FILES_PHYSICAL_CONSTANTS
  ${PROTOBUF_OUTPUT_FILES_BASE_PATH}/PhysicalConstants)

LIST(APPEND PROTOBUF_INPUT_FILES_PHYSICAL_CONSTANTS
  ${PROTOBUF_INPUT_FILES_BASE_PATH}/PhysicalConstants/Fundamental.proto)

# TODO: Not necessary for SET_SOURCE_FILES_PROPERTIES?
#ADD_SUBDIRECTORY(${PROTOBUF_INPUT_FILES_BASE_PATH}/PhysicalConstants
#  ${CMAKE_CURRENT_BINARY_DIR}/PhysicalConstants)

#SET_SOURCE_FILES_PROPERTIES(
#    ${PROTOBUF_INPUT_FILES_BASE_PATH}/PhysicalConstants
#    ${PROTOBUF_OUTPUT_FILES_PHYSICAL_CONSTANTS} PROPERTIES GENERATED TRUE)
 
# https://stackoverflow.com/questions/41430991/cmake-with-protobuf-file-in-subdirectory

FOREACH(proto ${PROTOBUF_INPUT_FILES_PHYSICAL_CONSTANTS})
  # TO_NATIVE_PATH mode converts cmake-style <path> into native path with
  # platform-specific slashes.
  FILE(TO_NATIVE_PATH ${proto} proto_native)

  IF(EXISTS ${proto_native})
    MESSAGE(${proto_native} " Exists")
  ELSE()
    MESSAGE(${proto_native} " Doesn't Exist")
  ENDIF()

  # cf. protoc --help protoc [OPTION] PROTO_FILES
  EXECUTE_PROCESS(COMMAND ${PROTOBUF_PROTOC_EXECUTABLE}
    # --proto_path=PATH Specify directory in which to search for imports. May be
    # specified multiple times; directories will be searched in order. If not
    # given, current working directory is used.
    # If not found in any of these directories, the --descriptor_set_in
    # descriptors will be checked for required proto file.
    --proto_path=${PROTOBUF_INPUT_FILES_BASE_PATH}/PhysicalConstants/
    # Generate C++ header and source
    --cpp_out=${PROTOBUF_OUTPUT_FILES_PHYSICAL_CONSTANTS} ${proto_native}
    RESULT_VARIABLE rv)

  # Optional, but that can show the user if something have gone wrong with the proto generation 
  IF(${rv})
    MESSAGE("Generation of data model returned ${rv} for proto ${proto_native}")
  ENDIF()
ENDFOREACH(proto)

ADD_LIBRARY(ProtoBinding
  ${PROTOBUF_OUTPUT_FILES_PHYSICAL_CONSTANTS}/Fundamental.pb.cc
  )

TARGET_LINK_LIBRARIES(ProtoBinding ${PROTOBUF_LIBRARIES})