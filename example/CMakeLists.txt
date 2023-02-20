###########################################################
##################### MEM EKF Example #####################
###########################################################
find_package(Python3 COMPONENTS Interpreter Development)
if (${Python3_FOUND})
  include_directories(${Python3_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "Python3 not found, please install it.")
endif()

add_executable(${EXAMPLE_PROJECT_NAME}
  main_example.cpp
  trajectory_generation.cpp
  matplotlibcpp.hpp
)

include_directories(${EXAMPLE_PROJECT_NAME}
  PUBLIC
  $<TARGET_PROPERTY:mem_ekf,INTERFACE_INCLUDE_DIRECTORIES>
  ${Python3_INCLUDE_DIRS}
)

target_link_libraries(${EXAMPLE_PROJECT_NAME}
  mem_ekf
  ${Python3_LIBRARIES}
)
