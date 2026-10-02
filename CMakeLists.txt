# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.16)
project(hydrocal VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(ENABLE_WARNINGS "Enable a common set of compiler warnings" ON)
option(WARNINGS_AS_ERRORS "Treat compiler warnings as errors (-Werror)" ON)
# FLINT is looked up through pkg-config, which resolves against the *host*
# even in a cross-compilation.  Turning this off keeps a host libflint from
# being linked into a foreign-target binary; the portable three-stage 1F1
# fallback is used instead.
option(USE_FLINT "Use FLINT for the exact 1F1/2F1 engines when available" ON)

# ---------------------------------------------------------------------------
# Warning flags (applied to all targets defined below)
# ---------------------------------------------------------------------------
if(ENABLE_WARNINGS)
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wshadow)
    if(WARNINGS_AS_ERRORS)
      add_compile_options(-Werror)
    endif()
  elseif(MSVC)
    add_compile_options(/W4)
    if(WARNINGS_AS_ERRORS)
      add_compile_options(/WX)
    endif()
  endif()
endif()

# ---------------------------------------------------------------------------
# Build-time generator for buildinfo.h
# ---------------------------------------------------------------------------
# Replacement for the old `svnversion -n` mechanism.
#
# HYDROCAL_REVISION continues the Subversion revision numbering, so numbers
# cited in the literature stay unique and ordered across the SVN -> git
# transition.  The counter lives in the REVISION file at the top of the
# source tree and is bumped by the pre-commit hook in .githooks/, which
# stages the new value together with the commit itself.  REVISION_BASE is
# only a fallback for checkouts where the file is missing (e.g. someone
# builds from a partial copy); in that case the value is reconstructed as
#   REVISION_BASE + (commits on the main line) - 1
# using `rev-list --first-parent`, so merges add one step rather than one per
# merged commit.
#
# HYDROCAL_GITID additionally records `git describe` output, which pins the
# exact commit (and shows "-dirty" for uncommitted changes) — useful when a
# bug report needs to name a precise build.
set(REVISION_BASE 2122 CACHE STRING
    "Fallback revision when the REVISION file is absent (last SVN revision)")

set(BUILDINFO_HEADER "${CMAKE_CURRENT_BINARY_DIR}/buildinfo.h")
file(WRITE "${CMAKE_BINARY_DIR}/gen_buildinfo.cmake" [=[
# Primary source: the counter file maintained by the pre-commit hook.
set(REV "")
set(COUNTER "${SRC_DIR}/REVISION")
if(EXISTS "${COUNTER}")
  file(READ "${COUNTER}" RAW)
  string(STRIP "${RAW}" RAW)
  if(RAW MATCHES "^[0-9]+$")
    set(REV "${RAW}")
  else()
    message(WARNING "REVISION does not contain a plain integer: '${RAW}'")
  endif()
endif()

# Fallback: reconstruct from the history when the counter file is missing.
if(REV STREQUAL "")
  execute_process(
    COMMAND git -C "${SRC_DIR}" rev-list --first-parent --count HEAD
    OUTPUT_VARIABLE NCOMMITS
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(NOT NCOMMITS MATCHES "^[0-9]+$")
    set(NCOMMITS 0)
  endif()
  math(EXPR REV "${REV_BASE} + ${NCOMMITS} - 1")
  if(REV LESS 0)
    set(REV 0)
  endif()
endif()

execute_process(
  COMMAND git -C "${SRC_DIR}" describe --tags --always --dirty
  OUTPUT_VARIABLE GITID
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
if(GITID STREQUAL "")
  set(GITID "unknown")
endif()
execute_process(
  COMMAND git -C "${SRC_DIR}" log -1 --format=%cd --date=short
  OUTPUT_VARIABLE REVDATE
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET
)
set(TEMP "${HEADER}.tmp")
file(WRITE "${TEMP}"
  "#ifndef HYDROCAL_BUILDINFO_H\n"
  "#define HYDROCAL_BUILDINFO_H\n"
  "#define HYDROCAL_VERSION \"${VERSION}\"\n"
  "#define HYDROCAL_REVISION \"${REV}\"\n"
  "#define HYDROCAL_GITID \"${GITID}\"\n"
  "#define HYDROCAL_REVISION_DATE \"${REVDATE}\"\n"
  "#endif\n")
if(EXISTS "${HEADER}")
  file(READ "${HEADER}" OLD)
  file(READ "${TEMP}" NEW)
  if(OLD STREQUAL NEW)
    file(REMOVE "${TEMP}")
    return()
  endif()
endif()
file(RENAME "${TEMP}" "${HEADER}")
]=])

add_custom_target(gen_buildinfo ALL
  COMMAND ${CMAKE_COMMAND} -D SRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
                          -D VERSION="${PROJECT_VERSION}"
                          -D REV_BASE="${REVISION_BASE}"
                          -D HEADER="${BUILDINFO_HEADER}"
                          -P "${CMAKE_BINARY_DIR}/gen_buildinfo.cmake"
  BYPRODUCTS "${BUILDINFO_HEADER}"
  COMMENT "Updating buildinfo.h"
)

# Git does not track hook configuration, so a fresh clone has no active
# hooks and the REVISION counter would silently stop advancing.  This target
# points git at the versioned .githooks/ directory.
find_program(GIT_EXECUTABLE git)
if(GIT_EXECUTABLE)
  add_custom_target(setup-hooks
    COMMAND ${GIT_EXECUTABLE} config core.hooksPath .githooks
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Activating .githooks (revision counter hook)"
  )
else()
  add_custom_target(setup-hooks
    COMMAND ${CMAKE_COMMAND} -E echo "git not found — run manually: git config core.hooksPath .githooks"
    COMMENT "git not found"
  )
endif()

# ---------------------------------------------------------------------------
# Library sources (all computational modules, no main)
# ---------------------------------------------------------------------------
set(HYDROCAL_SOURCES
  src/autoso1.cxx
  src/BetheBloch.cxx
  src/burgess.cxx
  src/clebsch.cxx
  src/Cooler.cxx
  src/CoolerDeconvolve.cxx
  src/CoolerFractions.cxx
  src/CoolerMonteCarlo.cxx
  src/DiracAngle.cxx
  src/DiracRate.cxx
  src/DiracTrans.cxx
  src/DiracRR.cxx
  src/DRmodelspe.cxx
  src/Faddeeva_w.cxx
  src/fele.cxx
  src/fieldion.cxx
  src/fieldion2.cxx
  src/FranckCondon.cxx
  src/gaussint.cxx
  src/hydromath.cxx
  src/hypergeometric.cxx
  src/kinema.cxx
  src/lifetime.cxx
  src/matrix.cxx
  src/numerov.cxx
  src/osci.cxx
  src/peakfunctions.cxx
  src/PIPE-MonteCarlo.cxx
  src/polari.cxx
  src/radrate.cxx
  src/readadas.cxx
  src/readstrz.cxx
  src/readxsec.cxx
  src/RRDRconvolve.cxx
  src/RRDRratecoef.cxx
  src/RRDRxsec.cxx
  src/RRphotons.cxx
  src/rydserie.cxx
  src/Stark.cxx
  src/storagetime.cxx
  src/stripping.cxx
  src/TOPbase.cxx
)

set(HYDROCAL_HEADERS
  src/autoso1.h
  src/BetheBloch.h
  src/burgess.h
  src/clebsch.h
  src/Cooler.h
  src/CoolerDeconvolve.h
  src/CoolerFractions.h
  src/CoolerMonteCarlo.h
  src/DiracAngle.h
  src/DiracRate.h
  src/DiracTrans.h
  src/DiracRR.h
  src/DRmodelspe.h
  src/Faddeeva_w.h
  src/fele.h
  src/fieldion.h
  src/fieldion2.h
  src/FranckCondon.h
  src/gaussint.h
  src/hydroconst.h
  src/hydromath.h
  src/hypergeometric.h
  src/IchiharaEichlerRRdata.h
  src/kinema.h
  src/lifetime.h
  src/matrix.h
  src/numerov.h
  src/osci.h
  src/peakfunctions.h
  src/PIPE-MonteCarlo.h
  src/polari.h
  src/radrate.h
  src/readadas.h
  src/readstrz.h
  src/readxsec.h
  src/RRDRconvolve.h
  src/RRDRratecoef.h
  src/RRDRxsec.h
  src/RRphotons.h
  src/rydserie.h
  src/Stark.h
  src/storagetime.h
  src/stripping.h
buildinfo.h
  src/TOPbase.h
)

# ---------------------------------------------------------------------------
# Library
# ---------------------------------------------------------------------------
add_library(hydrocal_lib STATIC ${HYDROCAL_SOURCES} ${HYDROCAL_HEADERS})
add_dependencies(hydrocal_lib gen_buildinfo)
target_include_directories(hydrocal_lib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/src" "${CMAKE_CURRENT_BINARY_DIR}")

# OpenMP
find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  target_link_libraries(hydrocal_lib PUBLIC OpenMP::OpenMP_CXX)
  target_compile_options(hydrocal_lib PRIVATE -fopenmp)
else()
  message(STATUS "OpenMP not found — compiling without parallel support")
endif()

# Boost (optional, header-only). Detection is done in the sources themselves via
# __has_include(<boost/...>); CMake's only job is to put Boost's headers on the
# include path so that self-detection can find them (needed for non-standard
# Boost install prefixes; system Boost in /usr/include is found without help).
find_package(Boost CONFIG QUIET)
if(Boost_FOUND)
  message(STATUS "Boost ${Boost_VERSION} found — Boost-based routines auto-enabled")
  target_link_libraries(hydrocal_lib PUBLIC Boost::headers)
else()
  message(STATUS "Boost not found — Boost-based routines disabled")
endif()

# FLINT (optional, non-header-only). Enables exact arbitrary-precision
# evaluation of the confluent hypergeometric functions 1F1 and 2F1 that
# underlie the closed-form (analytic) Dirac continuum wavefunction and
# pertaining matrix elements, respectively The sources
# self-detect the headers via __has_include just like Boost; CMake's job is
# to locate the library via pkg-config, link it, and advertise the
# HYDROCAL_HAVE_FLINT define.  Without FLINT the portable three-stage
# (double Taylor, quad Taylor, Poincare asymptotic) 1F1 fallback is used.
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND AND USE_FLINT)
  pkg_check_modules(FLINT QUIET IMPORTED_TARGET flint)
endif()
if(FLINT_FOUND)
  message(STATUS "FLINT ${FLINT_VERSION} found — exact 1F1 and 2F1 engines enabled")
  target_link_libraries(hydrocal_lib PUBLIC PkgConfig::FLINT)
  target_compile_definitions(hydrocal_lib PUBLIC HYDROCAL_HAVE_FLINT=1)
else()
  message(STATUS "FLINT not found — portable three-stage 1F1 fallback in use")
endif()

# Math library
target_link_libraries(hydrocal_lib PUBLIC m)

# ---------------------------------------------------------------------------
# Console executable
# ---------------------------------------------------------------------------
add_executable(hydrocal src/hydrocal.cxx)
target_link_libraries(hydrocal PUBLIC hydrocal_lib)

# ---------------------------------------------------------------------------
# GUI executable (requires FLTK + Cairo)
# ---------------------------------------------------------------------------
find_package(FLTK QUIET)
if(FLTK_FOUND)
  execute_process(COMMAND fltk-config --version OUTPUT_VARIABLE FLTK_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  message(STATUS "FLTK ${FLTK_VERSION} found")
else()
  message(STATUS "FLTK not found")
endif()

find_package(PkgConfig QUIET)
set(HAVE_CAIRO FALSE)
if(PkgConfig_FOUND)
  pkg_check_modules(CAIRO cairo QUIET IMPORTED_TARGET)
  if(CAIRO_FOUND)
    set(HAVE_CAIRO TRUE)
  endif()
endif()
if(HAVE_CAIRO)
  message(STATUS "Cairo found")
else()
  message(STATUS "Cairo not found")
endif()
if(FLTK_FOUND AND HAVE_CAIRO)
  add_executable(hydrocalGUI src/hydrocalGUI.cxx)
  target_link_libraries(hydrocalGUI PUBLIC hydrocal_lib ${FLTK_LIBRARIES})
  target_include_directories(hydrocalGUI PRIVATE ${FLTK_INCLUDE_DIRS} ${CAIRO_INCLUDE_DIRS})
  if(TARGET PkgConfig::CAIRO)
    target_link_libraries(hydrocalGUI PUBLIC PkgConfig::CAIRO)
  else()
    target_link_libraries(hydrocalGUI PUBLIC ${CAIRO_LIBRARIES})
    target_include_directories(hydrocalGUI PRIVATE ${CAIRO_INCLUDE_DIRS})
  endif()
  if(TARGET fltk::fltk)
    target_link_libraries(hydrocalGUI PUBLIC fltk::fltk)
  endif()
  message(STATUS "GUI enabled (FLTK + Cairo found)")
else()
  message(STATUS "GUI disabled (FLTK: ${FLTK_FOUND}, Cairo: ${HAVE_CAIRO})")
endif()

# ---------------------------------------------------------------------------
# Tests (CTest)
# ---------------------------------------------------------------------------
enable_testing()
add_test(
  NAME U92RR
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/U92RR.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/U92RR.result"
)
add_test(
  NAME Xe53lorentz
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/Xe53lorentz.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/Xe53lorentz.result"
)
add_test(
  NAME CG
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/CG.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/CG.result"
)
add_test(
  NAME ThreeJ
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/ThreeJ.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/ThreeJ.result"
)
add_test(
  NAME SixJ
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/SixJ.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/SixJ.result"
)
add_test(
  NAME NineJ
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/NineJ.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/NineJ.result"
)
add_test(
  NAME OsciBB
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/OsciBB.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/OsciBB.result"
)
add_test(
  NAME OsciBC
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/OsciBC.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/OsciBC.result"
)
add_test(
  NAME Hlifetimes
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/Hlifetimes.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/Hlifetimes.result"
)
add_test(
  NAME U91DiracRates
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/U91DiracRates.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/U91DiracRates.result"
)
add_test(
  NAME Pb82RRphotons
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/Pb82RRphotons.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/Pb82RRphotons.result"
)
set_tests_properties(Pb82RRphotons PROPERTIES TIMEOUT 300)
add_test(
  NAME Ds3p32DiracRR
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/Ds3p32DiracRR.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/Ds3p32DiracRR.result"
)
add_test(
  NAME S3soft
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/S3soft.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/S3soft.result"
)
add_test(
  NAME O2cooler
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/O2cooler.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/O2cooler.result"
)
add_test(
  NAME Ne2deconv
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/Ne2deconv.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/Ne2deconv.result"
)
add_test(
  NAME CNplus
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/run_test.sh
  "${CMAKE_CURRENT_SOURCE_DIR}/doc/examples/CNplus.hcin"
  "${CMAKE_CURRENT_BINARY_DIR}/CNplus.result"
)

# ---------------------------------------------------------------------------
# Symbolic links in ~/bin (original makefile behavior)
# ---------------------------------------------------------------------------
set(BIN_DIR "$ENV{HOME}/bin")
set(HYDROCAL_EXE "$<TARGET_FILE:hydrocal>")
install(CODE "file(MAKE_DIRECTORY \"${BIN_DIR}\")")
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
  \"${HYDROCAL_EXE}\" \"${BIN_DIR}/hydrocal\")")
if(TARGET hydrocalGUI)
  set(HYDROCALGUI_EXE "$<TARGET_FILE:hydrocalGUI>")
  install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    \"${HYDROCALGUI_EXE}\" \"${BIN_DIR}/hydrocalgui\")")
endif()

# ---------------------------------------------------------------------------
# Custom targets matching original makefile style
# ---------------------------------------------------------------------------
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/run_tests.sh" "#!/bin/sh
# Wrapper for ctest that supports: make run-tests TEST=<name>
if [ -n \"\$TEST\" ]; then
  ${CMAKE_CTEST_COMMAND} --output-on-failure -R \"\$TEST\"
else
  ${CMAKE_CTEST_COMMAND} --output-on-failure
fi
")
file(CHMOD "${CMAKE_CURRENT_BINARY_DIR}/run_tests.sh" PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)

add_custom_target(run-tests
  COMMAND sh "${CMAKE_CURRENT_BINARY_DIR}/run_tests.sh"
  DEPENDS hydrocal
)
# Small config snippet to inject the correct INCLUDE_PATH for buildinfo.h
set(DOXYGEN_EXTRA "${CMAKE_CURRENT_BINARY_DIR}/doxygen_extra.cfg")
file(WRITE "${DOXYGEN_EXTRA}" "INCLUDE_PATH += ${CMAKE_CURRENT_BINARY_DIR}\n")

find_program(CLANG_FORMAT clang-format)
file(GLOB FORMAT_SOURCES src/*.cxx src/*.h)
if(CLANG_FORMAT)
  add_custom_target(format
    COMMAND clang-format -i ${FORMAT_SOURCES}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Formatting sources with clang-format"
  )
else()
  add_custom_target(format
    COMMAND ${CMAKE_COMMAND} -E echo "clang-format not found — skipping"
    COMMENT "clang-format not found — install with: sudo apt install clang-format"
  )
endif()
add_custom_target(doc-html
  COMMAND doxygen "${CMAKE_CURRENT_SOURCE_DIR}/doc/hydrocal_html.dox"
                  "${DOXYGEN_EXTRA}"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  DEPENDS gen_buildinfo
)
add_custom_target(doc-pdf
  COMMAND doxygen "${CMAKE_CURRENT_SOURCE_DIR}/doc/hydrocal_pdf.dox"
                  "${DOXYGEN_EXTRA}"
  COMMAND ${CMAKE_MAKE_PROGRAM} -C "${CMAKE_CURRENT_SOURCE_DIR}/doc/latex"
  COMMAND cp "${CMAKE_CURRENT_SOURCE_DIR}/doc/latex/refman.pdf"
             "${CMAKE_CURRENT_SOURCE_DIR}/doc/hydrocal_manual.pdf"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
  DEPENDS gen_buildinfo
)
