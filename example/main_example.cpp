#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <string>

#include "measurement.hpp"
#include "mem_ekf.hpp"
#include "../components/helpers/helper_functions.hpp"

#include "trajectory_generation.hpp"

constexpr auto state_size = 4u;
constexpr auto measurement_size = 2u;


/*************************** Define motion model ***************************/
namespace eot {
  /* x, y, vx, vy */
  class ModelCv : public MemEkf<state_size> {
    public:
      explicit ModelCv(const MemEkfCalibrations<state_size> & calibrations) 
        : MemEkf<state_size>(calibrations) {
        // TODO
      }

    protected:
      void UpdateKinematic(const double time_delta) {
        // Update helpers
        SetTransitionMatrix(time_delta);

        // Update kinematic
        state_.kinematic_state.state = transition_matrix_ * state_.kinematic_state.state;
        state_.kinematic_state.covariance = transition_matrix_ * state_.kinematic_state.covariance * transition_matrix_.transpose() + time_delta * c_kinematic_;
      }
    
    private:
      void SetTransitionMatrix(const double time_delta) {
        //x
        transition_matrix_(0u, 0u) = 1.0;
        transition_matrix_(0u, 2u) = time_delta;
        // y
        transition_matrix_(1u, 2u) = 1.0;
        transition_matrix_(1u, 3u) = time_delta;
        // vx
        transition_matrix_(2u, 1u) = 1.0;
        // vy
        transition_matrix_(3u, 3u) = 1.0;
      }
      StateMatrix transition_matrix_ = StateMatrix::Zero();
  };
} //  namespace eot

/*************************** Main ***************************/
int main() {
  const auto gt = GetGroundTruth();

  std::default_random_engine generator;
  std::poisson_distribution<int> distribution(5.0);

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> rand(0.0, 1.0);

  /* Define tracker object */
  eot::MemEkfCalibrations<state_size> calibrations;
  calibrations.multiplicative_noise_diagonal = {0.25, 0.25};
  calibrations.process_noise_kinematic_diagonal = {100.0, 100.0, 1.0, 1.0};
  calibrations.process_noise_extent_diagonal = {0.05, 0.001, 0.001};
  calibrations.initial_state.kinematic_state.state << 100.0, 100.0, 10.0, -17.0;
  std::array<double, 4u> kin_cov = {900.0, 900.0, 16.0, 16.0};
  calibrations.initial_state.kinematic_state.covariance = eot::ConvertDiagonalToMatrix(kin_cov);
  calibrations.initial_state.extent_state.ellipse.alpha = std::numbers::pi_v<double> / 3.0;
  calibrations.initial_state.extent_state.ellipse.l1 = 200.0;
  calibrations.initial_state.extent_state.ellipse.l2 = 90.0;
  std::array<double, 3u> ext_cov = {0.2, 400.0, 400.0};
  calibrations.initial_state.extent_state.covariance = eot::ConvertDiagonalToMatrix(ext_cov);

  eot::ModelCv mem_ekf_cv_tracker(calibrations);

  /* Run */
  for (auto index = 0u; index < gt.time_steps; index++) {
    // Select detctions number in step
    auto detections_number = distribution(generator);
    while (detections_number == 0)
      detections_number = distribution(generator);

    std::cout << "Time step: " << std::to_string(index) << ", " << std::to_string(detections_number) << " Measurements\n";

    // Generate noisy measurement
    std::vector<eot::MeasurementWithCovariance<measurement_size>> measurements(detections_number);
    for (auto & measurement : measurements) {
      std::array<double, 2u> h = {-1.0 + 2 * rand(gen), -1.0 + 2 * rand(gen)};
      while (std::hypot(h.at(0), h.at(1)) > 1.0)
        h = {-1.0 + 2 * rand(gen), -1.0 + 2 * rand(gen)};

      measurement.value(0u) = gt.center.at(index).at(0u)
        + h.at(0) * gt.size.at(index).first * std::cos(gt.orientation.at(index))
        + h.at(1) * gt.size.at(index).second * std::cos(gt.orientation.at(index)); // + NOISE
      measurement.value(1u) = gt.center.at(index).at(1u)
        + h.at(0) * gt.size.at(index).first * std::sin(gt.orientation.at(index))
        + h.at(1) * gt.size.at(index).second * std::sin(gt.orientation.at(index)); // + NOISE

      measurement.covariance = Eigen::Matrix<double, measurement_size, measurement_size>::Zero();
      measurement.covariance(0u, 0u) = 200.0;
      measurement.covariance(1u, 1u) = 8.0;
    }

    // Run algo
    mem_ekf_cv_tracker.Run(static_cast<double>(index) * 10.0, measurements);
  }

  return EXIT_SUCCESS;
}