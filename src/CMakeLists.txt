set(HEADER_LIST
  "${SciCore_SOURCE_DIR}/include/SciCore/BasicMath.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/Cheb.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/Cheb2.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/ChebAdaptive.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/Definitions.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/DCT.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/DCT2.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/IDECheb.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/Integration.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/IntegrationWeights.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/Parallel.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/Random.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/Serialization.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/Utility.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/extern/CompilerVersions.h"
  "${SciCore_SOURCE_DIR}/include/SciCore/extern/FunctionTraits.h"
)

add_library(SciCore SHARED
    BasicMath.cpp
    Cheb.cpp
    Random.cpp
    Utility.cpp
    ${HEADER_LIST})
add_library(SciCore::SciCore ALIAS SciCore)

set_target_properties(SciCore PROPERTIES
    SOVERSION ${PROJECT_VERSION_MAJOR}
    VERSION ${PROJECT_VERSION}
)

include(GenerateExportHeader)
generate_export_header(SciCore
    EXPORT_FILE_NAME ${CMAKE_CURRENT_SOURCE_DIR}/../include/SciCore/SciCore_export.h
)

target_link_libraries(SciCore PRIVATE BLAS::BLAS LAPACK::LAPACK)
target_link_libraries(SciCore PUBLIC Eigen3::Eigen cereal::cereal Taskflow)

target_include_directories(SciCore PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(SciCore PUBLIC cxx_std_20)

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    # Every app/lib using SciCore needs to be build for the same architecture,
    # otherwise there can be segfaults using the Eigen library.
    # Therefore -march=native must be public.
    target_compile_options(SciCore PUBLIC
    $<$<CXX_COMPILER_ID:GNU>:-march=native>
    $<$<CXX_COMPILER_ID:Clang>:-march=native>
    )

    # These flags can be activatet only within the SciCore
    target_compile_options(SciCore PRIVATE
    $<$<CXX_COMPILER_ID:GNU>:-fno-math-errno -ffinite-math-only -fno-trapping-math -freciprocal-math -fno-signaling-nans -fcx-limited-range>
    $<$<CXX_COMPILER_ID:Clang>:-fno-math-errno -ffinite-math-only -fno-trapping-math -freciprocal-math -Wno-c99-extensions>
    )
endif()

# IDEs should put the headers in a nice place
source_group(
  TREE "${PROJECT_SOURCE_DIR}/include"
  PREFIX "Header Files"
  FILES ${HEADER_LIST})
