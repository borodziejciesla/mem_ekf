#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <map>
#include <numbers>
#include <random>
#include <set>
#include <string>

#include "matplotlibcpp.hpp"

#include "measurement.hpp"
#include "mem_ekf.hpp"
#include "../components/helpers/helper_functions.hpp"

#include "csv_reader.hpp"
#include "trajectory_generation.hpp"

namespace plt = matplotlibcpp;

constexpr auto state_size = 4u;
constexpr auto measurement_size = 2u;

//--------------------------------------------------------------------------//
//--- helper function convert timepoint to usable timestamp
template <typename TP>
time_t to_time_t(TP tp) {
  auto sctp = std::chrono::time_point_cast<std::chrono::system_clock::duration>(tp - TP::clock::now() + std::chrono::system_clock::now());
  return std::chrono::system_clock::to_time_t(sctp);
}

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
  /************************** Define tracker object **************************/
  eot::MemEkfCalibrations<state_size> calibrations;
  calibrations.multiplicative_noise_diagonal = {0.25, 0.25};
  calibrations.process_noise_kinematic_diagonal = {100.0, 100.0, 1.0, 1.0};
  calibrations.process_noise_extent_diagonal = {0.025, 0.00000001, 0.00000001};
  calibrations.initial_state.kinematic_state.state << 0.0, 0.0, 0.0, 0.0;
  std::array<double, 4u> kin_cov = {10.0, 10.0, 10.0, 10.0};
  calibrations.initial_state.kinematic_state.covariance = eot::ConvertDiagonalToMatrix(kin_cov);
  calibrations.initial_state.extent_state.ellipse.alpha = 0.0;
  calibrations.initial_state.extent_state.ellipse.l1 = 2.0;
  calibrations.initial_state.extent_state.ellipse.l2 = 1.0;
  std::array<double, 3u> ext_cov = {std::numbers::pi_v<double>, 1.0, 1.0};
  calibrations.initial_state.extent_state.covariance = eot::ConvertDiagonalToMatrix(ext_cov);

  eot::ModelCv mem_ekf_cv_tracker(calibrations);

  /************************** Run **************************/
  std::vector<eot::ObjectState<state_size>> output_objects;
  std::vector<std::vector<eot::MeasurementWithCovariance<measurement_size>>> detections;
  
  // Sort scans
  std::string radar_data_path("/home/maciek/Downloads/eot_simulation-20230227T213418Z-001/eot_simulation/lidar");

  std::set<std::filesystem::path> sorted_by_name;

  for (auto & entry : std::filesystem::directory_iterator(radar_data_path))
    sorted_by_name.insert(entry.path());

  // Run 
  double timestamp = 0.0;
  for (auto & scene_path : sorted_by_name) {
    //const auto scene_path = entry.path().c_str();
    std::cout << scene_path.c_str() << "\n";

    CsvReader reader(scene_path.c_str());
    const auto data = reader.GetData();

    std::vector<eot::MeasurementWithCovariance<measurement_size>> measurements;

    for (auto detection_index = 1u; detection_index < data.size(); detection_index++) {
      eot::MeasurementWithCovariance<measurement_size> measurement;

      measurement.value(0u) = std::stod(data.at(detection_index).at(0u));
      measurement.value(1u) = std::stod(data.at(detection_index).at(1u));

      measurement.covariance(0u, 0u) = 0.1;//std::stod(data.at(detection_index).at(1u));
      measurement.covariance(1u, 1u) = 0.1;//std::stod(data.at(detection_index).at(3u));

      measurements.push_back(measurement);
    }
    detections.push_back(measurements);

    std::cout << "Time step: " << std::to_string(timestamp) << ", " << std::to_string(measurements.size()) << " Measurements\n";

    // Run tracker
    mem_ekf_cv_tracker.Run(timestamp, measurements);

    // Get output
    const auto object = mem_ekf_cv_tracker.GetEstimatedState();
    output_objects.push_back(object);

    std::cout << "alpha = " << object.extent_state.ellipse.alpha << ", l1 = " << object.extent_state.ellipse.l1 << ", l2 = " << object.extent_state.ellipse.l2 << "\n";

    timestamp += 0.1;
  }

  /************************** Plot outputs **************************/
  plt::figure_size(1200, 780);
  plt::xlabel("X [m]");
  plt::ylabel("Y [m]");

  std::string reference_data_path("/home/maciek/Downloads/eot_simulation-20230227T213418Z-001/eot_simulation/reference");

  std::vector<double> x_traj;
  std::vector<double> y_traj;
  for (auto & scene_path : std::filesystem::directory_iterator(reference_data_path)) {
    CsvReader reader(scene_path.path().c_str());
    const auto data = reader.GetData();

    x_traj.push_back(std::stod(data.at(1u).at(0u)));
    y_traj.push_back(std::stod(data.at(1u).at(1u)));

    // eot::Ellipse ellipse = {std::stod(data.at(1u).at(2u)), 2.35, 0.9};
    // const auto [x_ellips, y_ellipse] = CreateEllipse(ellipse, std::make_pair(std::stod(data.at(1u).at(0u)), std::stod(data.at(1u).at(1u))));
    // plt::plot(x_ellips, y_ellipse, "k");
  }

  plt::plot(x_traj, y_traj, "^");

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
  // std::vector<double> x_ref;
  // std::vector<double> y_ref;
  // for (auto index = 0u; index < gt.time_steps; index = index + 3u) {
  //   x_ref.push_back(gt.center.at(index).at(0u));
  //   y_ref.push_back(gt.center.at(index).at(1u));

  //   eot::Ellipse ellipse = {gt.orientation.at(index), gt.size.at(index).first, gt.size.at(index).second};
  //   const auto [x_ellips, y_ellipse] = CreateEllipse(ellipse, std::make_pair(gt.center.at(index).at(0u), gt.center.at(index).at(1u)));
  //   plt::plot(x_ellips, y_ellipse, "k");
  // }
  // plt::plot(x_ref, y_ref, "k*");

  plt::axis("equal");
  plt::grid(true);
  plt::show();

  // Velocity
  // std::vector<double> vx_ref;
  // std::vector<double> vy_ref;
  std::vector<double> vx_obj;
  std::vector<double> vy_obj;
  std::vector<double> speed;
  std::vector<double> idx;
  for (auto index = 0u; index < output_objects.size(); index = index + 1u) {
    // vx_ref.push_back(gt.velocity.at(index).first);
    // vy_ref.push_back(gt.velocity.at(index).second);

    vx_obj.push_back(output_objects.at(index).kinematic_state.state(2u));
    vy_obj.push_back(output_objects.at(index).kinematic_state.state(3u));
    speed.push_back(std::hypot(output_objects.at(index).kinematic_state.state(2u), output_objects.at(index).kinematic_state.state(3u)));

    idx.push_back(index);
  }
  // plt::plot(idx, vx_ref, "b");
  // plt::plot(idx, vy_ref, "b:");
  plt::plot(idx, vx_obj, "r");
  plt::plot(idx, vy_obj, "r:");
  plt::plot(idx, speed, "r--");
  plt::show();

  return EXIT_SUCCESS;
}