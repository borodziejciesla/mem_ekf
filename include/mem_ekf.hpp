#ifndef MEM_EKF_INCLUDE_MEM_EKF_HPP_
#define MEM_EKF_INCLUDE_MEM_EKF_HPP_

#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "mem_ekf_calibrations.hpp"
#include "object_state.hpp"

namespace eot {
  template <size_t state_size, size_t measurement_size>
  class MemEkf {
    public:
      using StateVector = Eigen::Vector<double, state_size>;
      using StateCovariance = Eigen::Matrix<double, state_size, state_size>;
      using MeasurementVector = Eigen::Vector<double, measurement_size>;
      using MeasurementCovariance = Eigen::Matrix<double, measurement_size, measurement_size>;

    public:
      explicit MemEkf(const MemEkfCalibrations & calibrations)
        : calibrations_{calibrations} {
        /* Set constants */
        // f_ = [1 0 0 0; 0 0 0 1; 0 1 0 0]
        f_(0u, 0u) = 1.0;
        f_(1u, 3u) = 1.0;
        f_(2u, 1u) = 1.0;
        // f_tilde_ = [1 0 0 0; 0 0 0 1; 0 0 1 0];
        f_tilde_(0u, 0u) = 1.0;
        f_tilde_(1u, 3u) = 1.0;
        f_tilde_(2u, 2u) = 1.0;
      }
      
      virtual ~MemEkf(void) = default;

      void Run(const double timestamp, const std::vector<MeasurementVector> & measurements) {
        // Set time delta
        const auto time_delta = timestamp - prev_timestamp_;
        // Run algorithm
        RunUpdateStep(time_delta);
        RunCorrectionStep(measurements);
      }

      const ObjectState<state_size> & GetEstimatedState(void) const {
        return state_;
      }

    protected:
      virtual void UpdateKinematic(const double time_delta) = 0;

    private:
      void RunUpdateStep(const double time_delta) {
        /* Update kinematic */
        UpdateKinematic(time_delta);
        /* Update extent state */
        UpdateExtent();
      }

      void RunCorrectionStep(const std::vector<MeasurementVector> & measurements) {
        for (const auto & measurement : measurements)
          MakeOneDetectionCorrection(measurement);
      }

      void MakeOneDetectionCorrection(const MeasurementVector & measurement) {
        SetHelperVariables();

        // Calculate moments for the kinematic state update
        // Udpate kinematic estimate


        // Construct pseudo-measurement for the shape update
        // Calculate moments for the shape update
        // Update shape
      }

      void SetHelperVariables(void) {
        const auto alpha = state_.extent_state.ellipse.alpha;
        const auto l1 = state_.extent_state.ellipse.l1;
        const auto l2 = state_.extent_state.ellipse.l2;

        const auto c = std::cos(alpha);
        const auto s = std::sin(alpha);

        // Set s_
        s_(0u, 0u) = c * l1;
        s_(0u, 1u) = -s * l2;
        s_(1u, 0u) = s * l1;
        s_(1u, 1u) = c * l2;
        // Set s1_
        s1_(0u, 0u) = s_(0u, 0u);
        s1_(0u, 1u) = s_(0u, 1u);
        // Set s2_
        s2_(0u, 0u) = s_(1u, 0u);
        s2_(0u, 1u) = s_(1u, 1u);
        // Set j1_
        j1_(0u, 0u) = -l1 * s;
        j1_(0u, 1u) = c;
        j1_(1u, 0u) = -l2 * c;
        j1_(1u, 2u) = -s;
        // Set j2_
        j2_(0u, 0u) = l1 * c;
        j2_(0u, 1u) = s;
        j2_(1u, 0u) = -l2 * s;
        j2_(1u, 2u) = c;
        // Set c_i_
        c_i_ = s_ * c_h_ * s_.transpose();
        // Set c_ii_

      }

      void UpdateExtent(void) {
        //
      }

      MemEkfCalibrations calibrations_;

      double prev_timestamp_ = 0.0;

      Eigen::Matrix<double, 3u, 4u> f_ = Eigen::Matrix<double, 3u, 4u>::Zero();
      Eigen::Matrix<double, 3u, 4u> f_tilde_ = Eigen::Matrix<double, 3u, 4u>::Zero();
      Eigen::Matrix<double, 2u, 2u> s_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 1u, 2u> s1_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 1u, 2u> s2_ = Eigen::Matrix<double, 1u, 2u>::Zero();
      Eigen::Matrix<double, 2u, 3u> j1_ = Eigen::Matrix<double, 3u, 4u>::Zero();
      Eigen::Matrix<double, 2u, 3u> j2_ = Eigen::Matrix<double, 3u, 4u>::Zero();
      Eigen::Matrix<double, 2u, 2u> c_i_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 2u, 2u> c_ii_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 3u, 3u> m_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 4u, 2u> c_mz_ = Eigen::Matrix<double, 4u, 2u>::Zero();
      Eigen::Matrix<double, 2u, 2u> c_z_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 2u, state_size> h_ = Eigen::Matrix<double, 2u, state_size>::Zero();

      Eigen::Matrix<double, 2u, 2u> c_h_ = Eigen::Matrix<double, 2u, 2u>::Zero();

      // Estimated state
      ObjectState<state_size> state_;
  };
} //  namespace eot

#endif  //  MEM_EKF_INCLUDE_MEM_EKF_HPP_
