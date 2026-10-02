# SPDX-License-Identifier: MIT
#
# CMake toolchain for cross-compiling hydrocal to 64-bit Windows with
# mingw-w64.  Usage:
#
#   cmake -B build-win -S . \
#         -DCMAKE_TOOLCHAIN_FILE=cmake/toolchain-mingw-x86_64.cmake \
#         -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-win
#
# The result is build-win/hydrocal.exe, a statically linked binary that needs
# no additional DLLs.  See packaging/windows/make-win-zip.sh, which drives
# this toolchain and assembles the distributable archive.

set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_SYSTEM_PROCESSOR x86_64)

set(TOOLCHAIN_PREFIX x86_64-w64-mingw32)

set(CMAKE_C_COMPILER   ${TOOLCHAIN_PREFIX}-gcc)
set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)
set(CMAKE_RC_COMPILER  ${TOOLCHAIN_PREFIX}-windres)

# Only look for headers and libraries inside the mingw sysroot, but always use
# host programs (windres, the compiler driver itself).
set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX})
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# Link the runtime statically so the executable is self-contained and can be
# shipped without any companion DLLs.
# Link the runtime and OpenMP statically so the executable is self-contained
# and can be shipped without any companion DLLs.  Note that a plain -static is
# required and sufficient: -fopenmp plus -static resolves libgomp to its
# archive, whereas any -l<lib> named explicitly would pick the import library.
set(CMAKE_EXE_LINKER_FLAGS_INIT "-static -static-libgcc -static-libstdc++")

# The optional third-party dependencies (Boost, FLINT, FLTK, Cairo) are not
# typically available for the mingw target; hydrocal then falls back to its
# built-in implementations.  FLINT is looked up via pkg-config, which resolves
# against the host regardless of the toolchain, so it must be disabled
# explicitly — otherwise a host libflint would end up in the link line.
set(USE_FLINT OFF CACHE BOOL "" FORCE)

# CMake's FindOpenMP picks the libgomp import library when targeting Windows,
# which would leave hydrocal.exe depending on libgomp-1.dll.
set(USE_STATIC_OPENMP ON CACHE BOOL "" FORCE)

# Boost is detected with __has_include and is normally absent for the mingw
# target.  Code paths guarded by useBOOST then leave their arguments unused by
# design — the program reports at runtime that Boost is needed for those
# features — so warnings must not be fatal for this build.
set(WARNINGS_AS_ERRORS OFF CACHE BOOL "" FORCE)