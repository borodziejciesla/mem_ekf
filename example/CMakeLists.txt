###########################################################
##################### MEM EKF Example #####################
###########################################################
add_executable(${EXAMPLE_PROJECT_NAME}
  main_example.cpp
  trajectory_generation.cpp
)

include_directories(${EXAMPLE_PROJECT_NAME}
  PUBLIC
  $<TARGET_PROPERTY:mem_ekf,INTERFACE_INCLUDE_DIRECTORIES>
)

target_link_libraries(${EXAMPLE_PROJECT_NAME}
  mem_ekf
)
