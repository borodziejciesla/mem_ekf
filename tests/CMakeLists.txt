################################
# GTest
################################
project(googletest-git NONE)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        58d77fa  # release-1.10.0
)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

################################
# Tests
################################
add_executable(${TEST_PROJECT_NAME}
  mem_ekf_tests.cpp
)

target_link_libraries(${TEST_PROJECT_NAME}
  gtest 
  gtest_main
      
  mem_ekf
)
