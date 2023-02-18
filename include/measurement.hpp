#ifndef MEM_EKF_INCLUDE_MEASUREMENT_HPP_
#define MEM_EKF_INCLUDE_MEASUREMENT_HPP_

#include <Eigen/Dense>

namespace eot {
  template <size_t measurement_size>
  struct MeasurementWithCovariance {
    Eigen::Vector<double, measurement_size> value = Eigen::Vector<double, measurement_size>::Zero();
    Eigen::Matrix<double, measurement_size, measurement_size> covariance = Eigen::Matrix<double, measurement_size, measurement_size> ::Zero();
  };
} //  namespace eot

#endif  //  MEM_EKF_INCLUDE_MEASUREMENT_HPP_
