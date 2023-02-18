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

  template <size_t size>
  struct KinematicStateWithCovariance {
    Eigen::Vector<double, size> state;
    Eigen::Matrix<double, state, state> covariance = Eigen::Matrix<double, state, state>::Zero();
  };

  template <size_t size>
  struct ObjectState {
    EllipseWithCovariance extent_state;
    KinematicStateWithCovariance<size> kinematic_state;
  };
} //  namespace eot

#endif  //  MEM_EKF_INCLUDE_OBJECT_STATE_HPP_
