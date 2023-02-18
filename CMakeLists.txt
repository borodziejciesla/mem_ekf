cmake_minimum_required (VERSION 3.22)

project("mem_ekf" VERSION 0.0.1
  DESCRIPTION "C++ implementation of MEM-EKF tracker"
  HOMEPAGE_URL "https://github.com/borodziejciesla/mem_ekf")

set(TEST_PROJECT_NAME ${PROJECT_NAME}_tests)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_TESTS "Build tests" ON)


find_package(Eigen3)

add_subdirectory(components)

target_include_directories(${PROJECT_NAME}
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<TARGET_PROPERTY:eigen,INTERFACE_INCLUDE_DIRECTORIES>
)

# target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)

# locations are provided by GNUInstallDirs
# install(
#   TARGETS ${PROJECT_NAME}
#   EXPORT ${PROJECT_NAME}_Targets
#   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
# )

include(CMakePackageConfigHelpers)
write_basic_package_version_file("mem_ekfConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION
  ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
)

# install(
#   EXPORT ${PROJECT_NAME}_Targets
#   FILE ${PROJECT_NAME}Targets.cmake
#   NAMESPACE ${PROJECT_NAME}::
#   DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake
# )

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include DESTINATION include)

set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")

###############################################
################ Unit Tests ###################
###############################################
if(BUILD_TESTS)    
    project(${TEST_PROJECT_NAME})

    add_subdirectory(tests)    
    enable_testing()
    add_test(${TEST_PROJECT_NAME} ${TEST_PROJECT_NAME})    
endif() #   BUILD_TESTS
