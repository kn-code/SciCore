# Testing library
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest
  GIT_TAG v1.14.0
)

include(GoogleTest)

# If we just have set(INSTALL_GTEST OFF), then also Eigen is not installed (?)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_executable(BasicMathTest BasicMathTest.cpp)
target_link_libraries(BasicMathTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(BasicMathTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(BasicMathTest)

add_executable(ChebTest ChebTest.cpp)
target_link_libraries(ChebTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(ChebTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(ChebTest)

add_executable(Cheb2Test Cheb2Test.cpp)
target_link_libraries(Cheb2Test PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(Cheb2Test PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(Cheb2Test)

add_executable(ChebAdaptiveTest ChebAdaptiveTest.cpp)
target_link_libraries(ChebAdaptiveTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(ChebAdaptiveTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(ChebAdaptiveTest)

add_executable(DCTTest DCTTest.cpp)
target_link_libraries(DCTTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(DCTTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(DCTTest)

add_executable(IDEChebTest IDEChebTest.cpp)
target_link_libraries(IDEChebTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(IDEChebTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(IDEChebTest)

add_executable(IntegrationTest IntegrationTest.cpp)
target_link_libraries(IntegrationTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(IntegrationTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(IntegrationTest)

add_executable(RandomTest RandomTest.cpp)
target_link_libraries(RandomTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(RandomTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(RandomTest)

add_executable(SerializationTest SerializationTest.cpp)
target_link_libraries(SerializationTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(SerializationTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(SerializationTest)

add_executable(UtilityTest UtilityTest.cpp)
target_link_libraries(UtilityTest PRIVATE SciCore GTest::gtest_main)
target_compile_definitions(UtilityTest PRIVATE ${SCICORE_COMPILE_DEFINITIONS})
gtest_discover_tests(UtilityTest)

# Define custom target for building tests
add_custom_target(build_tests COMMENT "Build all the test executables.")
add_dependencies(build_tests BasicMathTest)
add_dependencies(build_tests ChebTest)
add_dependencies(build_tests Cheb2Test)
add_dependencies(build_tests ChebAdaptiveTest)
add_dependencies(build_tests DCTTest)
add_dependencies(build_tests IDEChebTest)
add_dependencies(build_tests IntegrationTest)
add_dependencies(build_tests RandomTest)
add_dependencies(build_tests SerializationTest)
add_dependencies(build_tests UtilityTest)