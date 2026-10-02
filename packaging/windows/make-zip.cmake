# SPDX-License-Identifier: MIT
#
# Assemble hydrocal_win.zip from an already configured and built tree.
# Invoked by the "win-zip" target:
#
#   cmake -DBUILD_DIR=... -DWINPROJ_DIR=... -DHYDROCAL_EXE=...
#         -DSOURCE_DIR=... -P packaging/windows/make-zip.cmake
#
# The archive is written to <source root>/hydrocal_win.zip, which is ignored
# by git on purpose: it is published as a GitHub Release asset, not committed.

foreach(_required BUILD_DIR WINPROJ_DIR HYDROCAL_EXE SOURCE_DIR)
  if(NOT DEFINED ${_required})
    message(FATAL_ERROR "make-zip.cmake: -D${_required} is required")
  endif()
endforeach()

set(_stage "${BUILD_DIR}/win-zip")
file(REMOVE_RECURSE "${_stage}")
file(MAKE_DIRECTORY "${_stage}")

# The console binary is shipped under the historic name hydrocal_win.exe.
get_filename_component(_exe_name "${HYDROCAL_EXE}" NAME)
file(COPY "${HYDROCAL_EXE}" DESTINATION "${_stage}")
file(RENAME "${_stage}/${_exe_name}" "${_stage}/hydrocal_win.exe")

set(_contents
    hydrocal_win.exe
    hydrocal.sln
    hydrocal.vcxproj
    ReadMe.txt
    buildinfo.h)

set(_paths "")
foreach(_item IN LISTS _contents)
  if(_item STREQUAL "hydrocal_win.exe")
    set(_src "${_stage}/${_item}")
  else()
    set(_src "${WINPROJ_DIR}/${_item}")
  endif()
  if(NOT EXISTS "${_src}")
    message(FATAL_ERROR
      "make-zip.cmake: missing ${_src} (build the 'win-project' target first)")
  endif()
  file(COPY "${_src}" DESTINATION "${_stage}")
  list(APPEND _paths "${_item}")
endforeach()

set(_archive "${SOURCE_DIR}/hydrocal_win.zip")
file(REMOVE "${_archive}")
# WORKING_DIRECTORY keeps the staging directory name out of the archive, so
# the files sit at the top level as users expect.
file(ARCHIVE_CREATE OUTPUT "${_archive}" FORMAT zip
     WORKING_DIRECTORY "${_stage}"
     PATHS ${_paths})

message(STATUS "Wrote ${_archive}")