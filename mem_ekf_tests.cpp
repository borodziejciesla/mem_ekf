#include "gtest/gtest.h"

#include "mem_ekf.hpp"

class ModelCv : public eot::MemEkf<4u, 2u> {
  protected:
    void UpdateKinematic(const double time_delta) {
      //
    }
}

/* Tests */
class MemEkfTests : public ::testing::Test
{
  protected:
    void SetUp(void) override {}
};

TEST_F(MemEkfTests, DummyTest)
{
    EXPECT_TRUE(true);
}
