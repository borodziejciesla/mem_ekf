#include "gtest/gtest.h"

#include "model_cv.hpp"

#include "mem_ekf.hpp"

/* Tests */
class MemEkfTests : public ::testing::Test
{
  protected:
    void SetUp(void) override {}
};

TEST_F(MemEkfTests, DummyTest)
{
  // Set calibrations
  eot::MemEkfCalibrations<4u> calibrations;
  calibrations.multiplicative_noise_diagonal = {0.25, 0.25};
  calibrations.process_noise_kinematic_diagonal = {100.0, 100.0, 1.0, 1.0};
  calibrations.process_noise_extent_diagonal = {0.05, 0.001, 0.001};

  // Create Filter object
  eot::ModelCv model_cv(calibrations);
  EXPECT_TRUE(true);
}
