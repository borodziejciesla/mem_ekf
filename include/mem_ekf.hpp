#ifndef MEM_EKF_INCLUDE_MEM_EKF_HPP_
#define MEM_EKF_INCLUDE_MEM_EKF_HPP_

#include <cmath>
#include <vector>

#include <Eigen/Dense>

#include "measurement.hpp"
#include "mem_ekf_calibrations.hpp"
#include "object_state.hpp"
#include "../components/helpers/helper_functions.hpp"

namespace eot {
  template <size_t state_size, size_t measurement_size = 2u>
  class MemEkf {
    public:
      using StateVector = Eigen::Vector<double, state_size>;
      using StateMatrix = Eigen::Matrix<double, state_size, state_size>;
      using Measurement = MeasurementWithCovariance<measurement_size>;
      using MeasurementVector = Eigen::Vector<double, measurement_size>;
      using MeasurementMatrix = Eigen::Matrix<double, measurement_size, measurement_size>;

    public:
      explicit MemEkf(const MemEkfCalibrations<state_size> & calibrations)
        : calibrations_{calibrations}
        , state_{calibrations_.initial_state}
        , c_kinematic_{ConvertDiagonalToMatrix(calibrations_.process_noise_kinematic_diagonal)}
        , c_h_{ConvertDiagonalToMatrix(calibrations_.multiplicative_noise_diagonal)}
        , c_extent_{ConvertDiagonalToMatrix(calibrations_.process_noise_extent_diagonal)} {
        /* Set constants */
        // f_ = [1 0 0 0; 0 0 0 1; 0 1 0 0]
        f_(0u, 0u) = 1.0;
        f_(1u, 3u) = 1.0;
        f_(2u, 1u) = 1.0;
        // f_tilde_ = [1 0 0 0; 0 0 0 1; 0 0 1 0];
        f_tilde_(0u, 0u) = 1.0;
        f_tilde_(1u, 3u) = 1.0;
        f_tilde_(2u, 2u) = 1.0;
        // h_ = [I(2x2), 0(2,state_size-2)]
        h_(0u, 0u) = 1.0;
        h_(1u, 1u) = 1.0;
      }
      
      virtual ~MemEkf(void) = default;

      void Run(const double timestamp, const std::vector<Measurement> & measurements) {
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
      
      ObjectState<state_size> state_;
      const Eigen::Matrix<double, state_size, state_size> c_kinematic_;

    private:
      void RunUpdateStep(const double time_delta) {
        /* Update kinematic */
        UpdateKinematic(time_delta);
        /* Update extent state */
        UpdateExtent();
      }

      void RunCorrectionStep(const std::vector<Measurement> & measurements) {
        for (const auto & measurement : measurements)
          MakeOneDetectionCorrection(measurement);
      }

      void MakeOneDetectionCorrection(const Measurement & measurement) {
        SetHelperVariables();
        const auto predicted_measurement = h_ * state_.kinematic_state.state;

        // Make corrections
        const auto cy = MakeKinematicCorrection(measurement, predicted_measurement);
        MakeExtentCorrection(measurement, predicted_measurement, cy);
      }

      Eigen::Matrix<double, 2u, 2u> MakeKinematicCorrection(const Measurement & measurement, const Eigen::Vector<double, 2u> & predicted_measurement) {
        // Calculate moments for the kinematic state update
        const auto cry = state_.kinematic_state.covariance * h_.transpose();
        const auto cy = h_ * state_.kinematic_state.covariance * h_.transpose() + c_i_ + c_ii_ + measurement.covariance;
        // Udpate kinematic estimate
        state_.kinematic_state.state += cry * cy.inverse() * (measurement.value - predicted_measurement);
        state_.kinematic_state.covariance -= cry * cy.inverse() * cry.transpose();
        // Force covariance symetry
        state_.kinematic_state.covariance = MakeMatrixSymetric<state_size>(state_.kinematic_state.covariance);

        return cy;
      }

      void MakeExtentCorrection(const Measurement & measurement, const Eigen::Vector<double, 2u> & predicted_measurement, const Eigen::Matrix<double, 2u, 2u> & cy) {
        // Construct pseudo-measurement for the shape update
        const Eigen::Vector2d innovation = measurement.value - predicted_measurement;
        Eigen::Vector4d innovation_prod;
        innovation_prod(0u) = std::pow(innovation(0u), 2);
        innovation_prod(1u) = innovation(0u) * innovation(1u);
        innovation_prod(2u) = innovation(0u) * innovation(1u);
        innovation_prod(3u) = std::pow(innovation(1u), 2);
        const Eigen::Vector3d yi = f_ * innovation_prod; 
        // Calculate moments for the shape update
        const Eigen::Vector3d yi_bar = f_ * cy.reshaped(4u, 1u);
        const Eigen::Matrix3d cp_y = state_.extent_state.covariance * m_.transpose();
        const Eigen::Matrix3d c_y = f_ * KroneckerProduct<2u, 2u, 2u, 2u>(cy, cy) * (f_ + f_tilde_).transpose();
        // Update shape
        const Eigen::Vector3d updated_ellipse_vector = ConvertEllipseToVector(state_.extent_state.ellipse) + cp_y * c_y.inverse() * (yi - yi_bar);
        state_.extent_state.ellipse = ConvertVectorToEllipse(updated_ellipse_vector);
        state_.extent_state.covariance -= static_cast<Eigen::Matrix3d>(cp_y * c_y.inverse()) * cp_y.transpose();
        // Force covariance symetry
        state_.extent_state.covariance = MakeMatrixSymetric<3u>(state_.extent_state.covariance);
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
        c_ii_(0u, 0u) = (state_.extent_state.covariance * j1_.transpose() * c_h_ * j1_).trace();
        c_ii_(0u, 1u) = (state_.extent_state.covariance * j2_.transpose() * c_h_ * j1_).trace();
        c_ii_(1u, 0u) = (state_.extent_state.covariance * j1_.transpose() * c_h_ * j2_).trace();
        c_ii_(1u, 1u) = (state_.extent_state.covariance * j2_.transpose() * c_h_ * j2_).trace();
        // set m_
        m_.row(0) = 2.0 * s1_ * c_h_ * j1_;
        m_.row(1) = 2.0 * s2_ * c_h_ * j2_;
        m_.row(2) = s1_ * c_h_ * j2_ + s2_ * c_h_ * j1_;
      }

      void UpdateExtent(void) {
        // elleipse = ellipse
        state_.extent_state.covariance += c_extent_;
      }

      MemEkfCalibrations<state_size> calibrations_;

      double prev_timestamp_ = 0.0;

      Eigen::Matrix<double, 3u, 4u> f_ = Eigen::Matrix<double, 3u, 4u>::Zero();
      Eigen::Matrix<double, 3u, 4u> f_tilde_ = Eigen::Matrix<double, 3u, 4u>::Zero();
      Eigen::Matrix<double, 2u, 2u> s_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 1u, 2u> s1_ = Eigen::Matrix<double, 1u, 2u>::Zero();
      Eigen::Matrix<double, 1u, 2u> s2_ = Eigen::Matrix<double, 1u, 2u>::Zero();
      Eigen::Matrix<double, 2u, 3u> j1_ = Eigen::Matrix<double, 2u, 3u>::Zero();
      Eigen::Matrix<double, 2u, 3u> j2_ = Eigen::Matrix<double, 2u, 3u>::Zero();
      Eigen::Matrix<double, 2u, 2u> c_i_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 2u, 2u> c_ii_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 3u, 3u> m_ = Eigen::Matrix<double, 3u, 3u>::Zero();
      Eigen::Matrix<double, 4u, 2u> c_mz_ = Eigen::Matrix<double, 4u, 2u>::Zero();
      Eigen::Matrix<double, 2u, 2u> c_z_ = Eigen::Matrix<double, 2u, 2u>::Zero();
      Eigen::Matrix<double, 2u, state_size> h_ = Eigen::Matrix<double, 2u, state_size>::Zero();

      const Eigen::Matrix<double, 2u, 2u> c_h_;
      const Eigen::Matrix3d c_extent_;
  };
} //  namespace eot

#endif  //  MEM_EKF_INCLUDE_MEM_EKF_HPP_
