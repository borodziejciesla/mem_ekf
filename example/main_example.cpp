#include <array>
#include <cmath>
#include <iostream>
#include <numbers>
#include <random>
#include <string>

#include "matplotlibcpp.hpp"

#include "measurement.hpp"
#include "mem_ekf.hpp"
#include "../components/helpers/helper_functions.hpp"

#include "trajectory_generation.hpp"

namespace plt = matplotlibcpp;

constexpr auto state_size = 4u;
constexpr auto measurement_size = 2u;

/*************************** Plot helpers ***************************/
std::pair<std::vector<double>, std::vector<double>> CreateEllipse(const eot::Ellipse & ellipse, const std::pair<double, double> & center) {
  std::vector<double> x;
  std::vector<double> y;

  for (double t = 0.0; t < 2.0 * std::numbers::pi_v<double>; t += 0.01) {
      double x_t = ellipse.l1 * std::cos(t) * std::cos(ellipse.alpha) - ellipse.l2 * std::sin(t) * std::sin(ellipse.alpha) + center.first;
      double y_t = ellipse.l1 * std::cos(t) * std::sin(ellipse.alpha) + ellipse.l2 * std::sin(t) * std::cos(ellipse.alpha) + center.second;
      x.push_back(x_t);
      y.push_back(y_t);
  }

  return std::make_pair(x, y);
}

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
        state_.kinematic_state.covariance = transition_matrix_ * state_.kinematic_state.covariance * transition_matrix_.transpose() + c_kinematic_;
      }
    
    private:
      void SetTransitionMatrix(const double time_delta) {
        //x
        transition_matrix_(0u, 0u) = 1.0;
        transition_matrix_(0u, 2u) = time_delta;
        // y
        transition_matrix_(1u, 1u) = 1.0;
        transition_matrix_(1u, 3u) = time_delta;
        // vx
        transition_matrix_(2u, 2u) = 1.0;
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

  std::normal_distribution<double> normal_distribution(0.0, 1.0);

  /************************** Define tracker object **************************/
  eot::MemEkfCalibrations<state_size> calibrations;
  calibrations.multiplicative_noise_diagonal = {0.25, 0.25};
  calibrations.process_noise_kinematic_diagonal = {100.0, 100.0, 1.0, 1.0};
  calibrations.process_noise_extent_diagonal = {0.05, 0.001, 0.001};
  calibrations.initial_state.kinematic_state.state << 100.0, 100.0, 10.0, -17.0;
  std::array<double, 4u> kin_cov = {900.0, 900.0, 16.0, 16.0};
  calibrations.initial_state.kinematic_state.covariance = eot::ConvertDiagonalToMatrix(kin_cov);
  calibrations.initial_state.extent_state.ellipse.alpha = -std::numbers::pi_v<double> / 3.0;
  calibrations.initial_state.extent_state.ellipse.l1 = 200.0;
  calibrations.initial_state.extent_state.ellipse.l2 = 90.0;
  std::array<double, 3u> ext_cov = {0.2, 400.0, 400.0};
  calibrations.initial_state.extent_state.covariance = eot::ConvertDiagonalToMatrix(ext_cov);

  eot::ModelCv mem_ekf_cv_tracker(calibrations);

  /************************** Run **************************/
  std::vector<eot::ObjectState<state_size>> output_objects;
  std::vector<std::vector<eot::MeasurementWithCovariance<measurement_size>>> detections;

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
        - h.at(1) * gt.size.at(index).second * std::sin(gt.orientation.at(index));// + 10.0 * normal_distribution(generator);
      measurement.value(1u) = gt.center.at(index).at(1u)
        + h.at(0) * gt.size.at(index).first * std::sin(gt.orientation.at(index))
        + h.at(1) * gt.size.at(index).second * std::cos(gt.orientation.at(index));// + 10.0 * normal_distribution(generator);

      measurement.covariance = Eigen::Matrix<double, measurement_size, measurement_size>::Zero();
      measurement.covariance(0u, 0u) = 200.0;
      measurement.covariance(1u, 1u) = 8.0;
    }

    detections.push_back(measurements);

    // Run algo
    mem_ekf_cv_tracker.Run(static_cast<double>(index) * 10.0, measurements);
    output_objects.push_back(mem_ekf_cv_tracker.GetEstimatedState());


    const auto object = mem_ekf_cv_tracker.GetEstimatedState();
    std::cout << "alpha = " << object.extent_state.ellipse.alpha << ", l1 = " << object.extent_state.ellipse.l1 << ", l2 = " << object.extent_state.ellipse.l2 << "\n";
    std::cout << "Center Error [m]: " << std::hypot(object.kinematic_state.state(0) - gt.center.at(index).at(0), object.kinematic_state.state(1) - gt.center.at(index).at(1)) << "\n";
  }

  /************************** Plot outputs **************************/
  plt::figure_size(1200, 780);

  plt::xlabel("X [m]");
  plt::ylabel("Y [m]");

  // Trajectory
  std::vector<double> x_traj;
  std::vector<double> y_traj;
  for (const auto & point : gt.center) {
    x_traj.push_back(point.at(0u));
    y_traj.push_back(point.at(1u));
  }

  std::map<std::string, std::string> keywords_traj;
  keywords_traj.insert(std::pair<std::string, std::string>("label", "Trajectory") );

  plt::plot(x_traj, y_traj, keywords_traj);

  // Detections
  std::vector<double> x_detections;
  std::vector<double> y_detections;
  for (auto index = 0u; index < detections.size(); index = index + 3u) {
    for (const auto & detection : detections.at(index)) {
      x_detections.push_back(detection.value(0u));
      y_detections.push_back(detection.value(1u));
    }
  }
  plt::scatter(x_detections, y_detections);

  // Objects Center
  std::vector<double> x_objects;
  std::vector<double> y_objects;
  for (auto index = 0u; index < output_objects.size(); index = index + 3u) {
    x_objects.push_back(output_objects.at(index).kinematic_state.state(0u));
    y_objects.push_back(output_objects.at(index).kinematic_state.state(1u));

    const auto [x_ellips, y_ellipse] = CreateEllipse(output_objects.at(index).extent_state.ellipse, std::make_pair(output_objects.at(index).kinematic_state.state(0u), output_objects.at(index).kinematic_state.state(1u)));
    plt::plot(x_ellips, y_ellipse, "r");
  }
  plt::plot(x_objects, y_objects, "r*");

  // Reference
  std::vector<double> x_ref;
  std::vector<double> y_ref;
  for (auto index = 0u; index < gt.time_steps; index = index + 3u) {
    x_ref.push_back(gt.center.at(index).at(0u));
    y_ref.push_back(gt.center.at(index).at(1u));

    eot::Ellipse ellipse = {gt.orientation.at(index), gt.size.at(index).first, gt.size.at(index).second};
    const auto [x_ellips, y_ellipse] = CreateEllipse(ellipse, std::make_pair(gt.center.at(index).at(0u), gt.center.at(index).at(1u)));
    plt::plot(x_ellips, y_ellipse, "k");
  }
  plt::plot(x_ref, y_ref, "k*");
  plt::show();

  // Velocity
  std::vector<double> vx_ref;
  std::vector<double> vy_ref;
  std::vector<double> vx_obj;
  std::vector<double> vy_obj;
  std::vector<double> idx;
  for (auto index = 0u; index < gt.time_steps; index = index + 1u) {
    vx_ref.push_back(gt.velocity.at(index).first);
    vy_ref.push_back(gt.velocity.at(index).second);

    vx_obj.push_back(output_objects.at(index).kinematic_state.state(2u));
    vy_obj.push_back(output_objects.at(index).kinematic_state.state(3u));

    idx.push_back(index);
  }
  plt::plot(idx, vx_ref, "b");
  plt::plot(idx, vy_ref, "b:");
  plt::plot(idx, vx_obj, "r");
  plt::plot(idx, vy_obj, "r:");
  plt::show();

  return EXIT_SUCCESS;
}