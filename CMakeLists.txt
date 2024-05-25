cmake_minimum_required(VERSION 3.25)

project(
    SciCore
    VERSION 0.9.0
    DESCRIPTION "Computational Science in C++"
    HOMEPAGE_URL "https://github.com/kn-code/SciCore"
    LANGUAGES CXX
)

# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

include(FetchContent)

set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

    set(CMAKE_CXX_STANDARD 20)
    set(CMAKE_CXX_EXTENSIONS OFF)
    set(CMAKE_COLOR_DIAGNOSTICS ON)

    set(CMAKE_POSITION_INDEPENDENT_CODE ON)
    set(CMAKE_CXX_VISIBILITY_PRESET hidden)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

    if(NOT APPLE)
        set(CMAKE_INSTALL_RPATH $ORIGIN)
    endif()

    # Useful for debugging/benchmarking are also -g -fno-omit-frame-pointer
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wnon-virtual-dtor -pedantic -Wcast-align -Wunused -Woverloaded-virtual -Wnull-dereference -Wdouble-promotion")
    endif ()

    set_property(GLOBAL PROPERTY USE_FOLDERS ON)

    set(FETCHCONTENT_QUIET OFF)

    include(CTest)
endif()

#
# Dependencies
#

# BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

#
# Eigen
#

find_package(Eigen3 3.4.0 QUIET)

# If Eigen was not found use FetchContent to download it
if(NOT Eigen3_FOUND)
    message(STATUS "Eigen3 not found. Fetching from Git repository...")

    FetchContent_Declare(
        Eigen3
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG        3.4.0
    )

    block (SCOPE_FOR VARIABLES)
        set(BUILD_TESTING OFF)
        FetchContent_MakeAvailable(Eigen3)
    endblock ()
else()
    # Print the found Eigen include directory
    message(STATUS "Eigen3 found: ${EIGEN3_INCLUDE_DIR}")
endif()

#
# Cereal
#

FetchContent_Declare(
    cereal
    GIT_REPOSITORY https://github.com/USCiLab/cereal.git
    GIT_TAG        v1.3.2
)

block (SCOPE_FOR VARIABLES)
    set(CEREAL_INSTALL ON)
    set(BUILD_DOC OFF)
    set(BUILD_SANDBOX OFF)
    set(SKIP_PERFORMANCE_COMPARISON ON)
    set(BUILD_TESTING OFF)
    FetchContent_MakeAvailable(cereal)
endblock ()

#
# Taskflow
#

FetchContent_Declare(
    Taskflow
    GIT_REPOSITORY https://github.com/taskflow/taskflow.git
    GIT_TAG        v3.7.0
)

block (SCOPE_FOR VARIABLES)
    set(TF_BUILD_TESTS OFF)
    set(TF_BUILD_EXAMPLES OFF)
    FetchContent_MakeAvailable(Taskflow)
endblock ()

#
# Build library and tests
#

add_subdirectory(src)

# Testing only available if this is the main app and BUILD_TESTING is defined
if((CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME) AND BUILD_TESTING)
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()

#
# Installation instructions
#

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

set(SciCore_INSTALL_CMAKEDIR
    ${CMAKE_INSTALL_LIBDIR}/cmake/SciCore
    CACHE STRING "Path to SciCore CMake files"
)

install(
    DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    TARGETS SciCore cereal
    EXPORT SciCoreTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        COMPONENT   SciCore_Runtime
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT   SciCore_Runtime
        NAMELINK_COMPONENT SciCore_Development
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        COMPONENT SciCore_Development
)

install(EXPORT SciCoreTargets
    DESTINATION ${SciCore_INSTALL_CMAKEDIR}
    NAMESPACE SciCore::
    FILE SciCoreTargets.cmake
    COMPONENT SciCore_Development
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SciCoreConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/Config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/SciCoreConfig.cmake
    INSTALL_DESTINATION ${SciCore_INSTALL_CMAKEDIR}
)

install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/SciCoreConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/SciCoreConfigVersion.cmake
    DESTINATION ${SciCore_INSTALL_CMAKEDIR}
)

install(
    DIRECTORY ${cereal_SOURCE_DIR}/include/cereal
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

export(EXPORT SciCoreTargets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/SciCoreTargets.cmake"
       NAMESPACE SciCore::)
