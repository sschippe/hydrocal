# SPDX-License-Identifier: MIT
#
# Fill in packaging/windows/ReadMe.txt.in once buildinfo.h is available.
# Invoked by the "win-project" target:
#
#   cmake -DBUILD_DIR=... -DWINPROJ_DIR=... -DSOURCE_DIR=...
#         -P packaging/windows/make-readme.cmake
#
# The version macros only exist inside the generated buildinfo.h, and that
# header is produced during the build, so this cannot be done at configure time.

foreach(_required BUILD_DIR WINPROJ_DIR SOURCE_DIR)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "make-readme.cmake: -D${_required} is required")
  endif()
endforeach()

set(_header "${BUILD_DIR}/buildinfo.h")
if(NOT EXISTS "${_header}")
  message(FATAL_ERROR "make-readme.cmake: ${_header} not found (build first)")
endif()
file(READ "${_header}" _buildinfo)

foreach(_macro VERSION REVISION GITID REVISION_DATE)
  if(_buildinfo MATCHES "#define[ \t]+HYDROCAL_${_macro}[ \t]+\"([^\"]*)\"")
    set(HYDROCAL_${_macro} "${CMAKE_MATCH_1}")
  else()
    message(FATAL_ERROR "make-readme.cmake: HYDROCAL_${_macro} not found in ${_header}")
  endif()
endforeach()

configure_file(
  "${SOURCE_DIR}/packaging/windows/ReadMe.txt.in"
  "${WINPROJ_DIR}/ReadMe.txt"
  @ONLY)

message(STATUS
  "Wrote ${WINPROJ_DIR}/ReadMe.txt "
  "(version ${HYDROCAL_VERSION}, revision ${HYDROCAL_REVISION})")