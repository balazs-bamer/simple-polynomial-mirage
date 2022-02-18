#include "Angle2apparentMirrorDepth.h"
#include "gtest/gtest.h"


constexpr float cgEpsilon = 0.0001f;

bool eq(double const aF1, double const aF2, double const aEpsilon) {
  return ::abs(aF1 - aF2) < aEpsilon;
}

bool eq(double const aF1, double const aF2) {
  return eq(aF1, aF2, cgEpsilon);
}

TEST(polynomApprox, x2few) {
  PolynomApprox poly;
  poly.init(std::vector({1.0, 4.0, 9.0}), std::vector({1.0, 2.0, 3.0}), 2u);
  EXPECT_TRUE(eq(poly.eval(0.0),  0.0));
  EXPECT_TRUE(eq(poly.eval(1.0),  0.0));
  EXPECT_TRUE(eq(poly.eval(4.0),  0.0));
 // EXPECT_TRUE(eq(poly.getRrmsError(), std::numeric_limits<double>::max()));
}

TEST(polynomApprox, x20) {
  PolynomApprox poly;
  poly.init(std::vector({0.0, 0.0, 0.0}), std::vector({1.0, 2.0, 3.0}), 2u);
  EXPECT_TRUE(eq(poly.eval(0.0),  0.0));
  EXPECT_TRUE(eq(poly.eval(4.0),  0.0));
  EXPECT_TRUE(eq(poly.eval(5.0),  0.0));
 // EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_21_27_30) {
  PolynomApprox poly;
  poly.init(std::vector({21.0, 27.0, 30.0}), std::vector({32.0, 48.0, 63.0}), 2u);
  EXPECT_TRUE(eq(poly.eval(32.0),  21.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(48.0),  27.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(63.0),  30.0, 0.1));
 // EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2) {
  PolynomApprox poly;
  poly.init(std::vector({1.0, 4.0, 9.0}), std::vector({1.0, 2.0, 3.0}), 2u);
  EXPECT_TRUE(eq(poly.eval(-1.0), 1.0));
  EXPECT_TRUE(eq(poly.eval(0.0),  0.0));
  EXPECT_TRUE(eq(poly.eval(4.0), 16.0));
  EXPECT_TRUE(eq(poly.eval(5.0), 25.0));
  EXPECT_TRUE(eq(poly.eval(6.0), 36.0));
 // EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_5) {
  PolynomApprox poly;
  poly.init(std::vector({6.0, 9.0, 14.0}), std::vector({1.0, 2.0, 3.0}), 2u);
  EXPECT_TRUE(eq(poly.eval(-1.0), 6.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(0.0),  5.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(4.0), 21.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(5.0), 30.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(6.0), 41.0, 0.1));
 // EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2off) {
  PolynomApprox poly;
  poly.init(std::vector({1.0, 4.0, 9.0, 15.9}), std::vector({1.0, 2.0, 3.0, 4.0}), 2u);
  EXPECT_TRUE(eq(poly.eval(0.0),  0.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(4.0), 16.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(5.0), 25.0, 0.3));
 // EXPECT_TRUE(eq(poly.getRrmsError(), 0.00059f));
}

TEST(polynomApprox, x2veryoff) {
  PolynomApprox poly;
  poly.init(std::vector({1.0, 3.0, 13.0, 13.0}), std::vector({1.0, 2.0, 3.0, 4.0}), 2u);
  EXPECT_FALSE(eq(poly.eval(0.0),  0.0));
  EXPECT_FALSE(eq(poly.eval(4.0), 16.0));
  EXPECT_FALSE(eq(poly.eval(5.0), 25.0));
 // EXPECT_TRUE(poly.getRrmsError() > 0.1f);
}

TEST(polynomApprox, 2x3_3x2_4x_5) {
  PolynomApprox poly;
  poly.init(std::vector({5.0, 14.0, 41.0, 98.0}), std::vector({0.0, 1.0, 2.0, 3.0}), 3u);
  EXPECT_TRUE(eq(poly.eval(0.5),  8.0));
  EXPECT_TRUE(eq(poly.eval(1.5), 24.5));
  EXPECT_TRUE(eq(poly.eval(2.5), 65.0));
 // EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(angle2apparentMirrorDepth, height2temp) {
  Angle2apparentMirrorDepth mirage(28);
  EXPECT_TRUE(eq(1.0, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(1.0))));
  EXPECT_TRUE(eq(0.5, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.5))));
  EXPECT_TRUE(eq(0.1, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.1))));
  EXPECT_TRUE(eq(0.05, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.05))));
  EXPECT_TRUE(eq(0.01, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.01))));
  EXPECT_TRUE(eq(0.005, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.005))));
}

TEST(angle2apparentMirrorDepth, temp2height) {
  Angle2apparentMirrorDepth mirage(28);
  EXPECT_TRUE(eq(0.0, mirage.getTempRiseAtHeight(mirage.getHeightAtTempRise(0.0))));
  EXPECT_TRUE(eq(8.0, mirage.getTempRiseAtHeight(mirage.getHeightAtTempRise(8.0))));
  EXPECT_TRUE(eq(16.0, mirage.getTempRiseAtHeight(mirage.getHeightAtTempRise(16.0))));
  EXPECT_TRUE(eq(24.0, mirage.getTempRiseAtHeight(mirage.getHeightAtTempRise(24.0))));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
