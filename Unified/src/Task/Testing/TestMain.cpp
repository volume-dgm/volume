#include "../Task/Task.h"
#include <cstdlib>
#include "../../../3rdparty/googletest/googletest/include/gtest/gtest.h"

#define SPACE_FROM_SETTINGS
const unsigned int defaultPolynomialOrder = 1;

#include "TestRun.inl"
#include "TestRunDecomposed.inl"

TEST(InitTest, TestRun) {
  int dimCount = 2;
  const char* settingsName = "test_configs/test1.xml";

  Task<Space2, defaultPolynomialOrder> task(settingsName);
  task.TestRun();
}

int main(int argc, char* argv[]) {
  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}