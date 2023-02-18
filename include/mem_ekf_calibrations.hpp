#ifndef MEM_EKF_INCLUDE_MEM_EKF_CALIBRATIONS_HPP_
#define MEM_EKF_INCLUDE_MEM_EKF_CALIBRATIONS_HPP_

#include <array>
namespace eot {
  template <size_t state_size>
  struct MemEkfCalibrations {
    std::array<double, 2u> multiplicative_noise_diagonal;             // Covariance of the multiplicative noise
    std::array<double, state_size> process_noise_kinematic_diagonal;  // Covariance of the process noise for the kinematic state
    std::array<double, 3u> process_noise_extent_diagonal;             // Covariance of the process noise for the shape parameters
  };
} //  namespace eot

#endif  //  MEM_EKF_INCLUDE_MEM_EKF_CALIBRATIONS_HPP_
