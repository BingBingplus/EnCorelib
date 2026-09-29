cmake_minimum_required(VERSION 3.16)
project(EnCorelib CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ────────────────────────────────────────────────────────────────────────────
#  CAPD dependency
#  Try system installation first; fall back to FetchContent from GitHub.
# ────────────────────────────────────────────────────────────────────────────
option(USE_SYSTEM_CAPD "Use system-installed CAPD (via capd-config)" ON)
option(USE_SYMENGINE   "Enable SymEngine for algebraic expression expansion" OFF)

if(USE_SYSTEM_CAPD)
    # Try pkg-config / capd-config
    find_program(CAPD_CONFIG NAMES capd-config)
    if(CAPD_CONFIG)
        execute_process(COMMAND ${CAPD_CONFIG} --cflags
                        OUTPUT_VARIABLE CAPD_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(COMMAND ${CAPD_CONFIG} --libs
                        OUTPUT_VARIABLE CAPD_LIBS   OUTPUT_STRIP_TRAILING_WHITESPACE)
        message(STATUS "Found CAPD via capd-config: ${CAPD_CONFIG}")
        add_library(capd_import INTERFACE)
        separate_arguments(CAPD_CFLAGS_LIST NATIVE_COMMAND "${CAPD_CFLAGS}")
        separate_arguments(CAPD_LIBS_LIST   NATIVE_COMMAND "${CAPD_LIBS}")
        target_compile_options(capd_import INTERFACE ${CAPD_CFLAGS_LIST})
        target_link_libraries( capd_import INTERFACE ${CAPD_LIBS_LIST})
        set(CAPD_FOUND TRUE)
    else()
        set(CAPD_FOUND FALSE)
        message(WARNING "capd-config not found; will try FetchContent.")
    endif()
else()
    set(CAPD_FOUND FALSE)
endif()

if(NOT CAPD_FOUND)
    # Download CAPD from GitHub
    include(FetchContent)
    message(STATUS "Fetching CAPD from GitHub (first build may take several minutes)...")
    FetchContent_Declare(
        capd
        GIT_REPOSITORY https://github.com/CAPDGroup/CAPD.git
        GIT_TAG        master
        GIT_SHALLOW    TRUE
    )
    # CAPD options: disable tests/examples to speed up
    set(CAPD_BUILD_TESTS     OFF CACHE BOOL "" FORCE)
    set(CAPD_BUILD_EXAMPLES  OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(capd)
    # CAPD CMake exports a target; name may vary – try capd or CAPD::capd
    if(TARGET capd)
        add_library(capd_import INTERFACE)
        target_link_libraries(capd_import INTERFACE capd)
        set(CAPD_FOUND TRUE)
    elseif(TARGET CAPD::capd)
        add_library(capd_import INTERFACE)
        target_link_libraries(capd_import INTERFACE CAPD::capd)
        set(CAPD_FOUND TRUE)
    else()
        message(FATAL_ERROR "Could not find a CAPD target after FetchContent. "
                            "Please install CAPD manually and set USE_SYSTEM_CAPD=ON.")
    endif()
endif()

# ────────────────────────────────────────────────────────────────────────────
#  Optional: SymEngine
# ────────────────────────────────────────────────────────────────────────────
if(USE_SYMENGINE)
    find_package(SymEngine REQUIRED)
    message(STATUS "SymEngine found: ${SYMENGINE_INCLUDE_DIRS}")
endif()

# ────────────────────────────────────────────────────────────────────────────
#  Build target
# ────────────────────────────────────────────────────────────────────────────
# -frounding-math is REQUIRED: without it the compiler may assume the
# default rounding mode, which breaks the directed rounding that every
# interval bound in this program relies on.
function(EnCorelib_configure tgt)
    target_link_libraries(${tgt} PRIVATE capd_import)
    target_compile_options(${tgt} PRIVATE -Wall -O2 -frounding-math)
    if(USE_SYMENGINE)
        target_compile_definitions(${tgt} PRIVATE HAVE_SYMENGINE)
        target_include_directories(${tgt} PRIVATE ${SYMENGINE_INCLUDE_DIRS})
        target_link_libraries(${tgt} PRIVATE symengine gmp mpfr mpc)
    endif()
endfunction()

add_executable(EnCorelib EnCorelib.cpp)
EnCorelib_configure(EnCorelib)

# Validated trajectory printer used for Figure 1 (top row).
add_executable(Encoretraj tools/trajectory.cpp)
EnCorelib_configure(Encoretraj)

# CAPD (Cr-Lohner, r = 3) reference used for the CAPD rows of Table 2.
add_executable(capdref repro/competitors/capd_compare.cpp)
EnCorelib_configure(capdref)

# ────────────────────────────────────────────────────────────────────────────
#  Self-test:  ctest    (or: bash repro/check.sh)
# ────────────────────────────────────────────────────────────────────────────
enable_testing()
add_test(NAME selftest
         COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/repro/check.sh
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

# ────────────────────────────────────────────────────────────────────────────
#  Install
# ────────────────────────────────────────────────────────────────────────────
install(TARGETS EnCorelib Encoretraj capdref DESTINATION bin)
install(FILES
    types.h symparse.h stepAB.h extend.h refine.h endcover.h boundary.h
    DESTINATION include/EnCorelib)
install(DIRECTORY examples/ DESTINATION share/EnCorelib/examples)
