#ifndef MEM_EKF_INCLUDE_OBJECT_STATE_HPP_
#define MEM_EKF_INCLUDE_OBJECT_STATE_HPP_

#include <Eigen/Dense>

#include "ellipse.hpp"

namespace eot {
  using EllipseCovariance = Eigen::Matrix<double, 3u, 3u>;

  struct EllipseWithCovariance {
    Ellipse ellipse;
    EllipseCovariance covariance = EllipseCovariance::Zero();
  };

  template <size_t state_size>
  struct KinematicStateWithCovariance {
    Eigen::Vector<double, state_size> state = Eigen::Vector<double, state_size>::Zero();
    Eigen::Matrix<double, state_size, state_size> covariance = Eigen::Matrix<double, state_size, state_size>::Zero();
  };

  template <size_t state_size>
  struct ObjectState {
    EllipseWithCovariance extent_state;
    KinematicStateWithCovariance<state_size> kinematic_state;
  };
} //  namespace eot

#endif  //  MEM_EKF_INCLUDE_OBJECT_STATE_HPP_
