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
  eot::MemEkfCalibrations<4u> calibrations;
  eot::ModelCv model_cv(calibrations);
  EXPECT_TRUE(true);
}
