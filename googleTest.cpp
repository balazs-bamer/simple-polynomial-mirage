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
  PolynomApprox poly({1.0f, 2.0f}, {1.0f, 4.0f, 9.0f}, 2u);
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(1.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f),  0.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), std::numeric_limits<double>::max()));
}

TEST(polynomApprox, x20) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f}, {0.0f, 0.0f, 0.0f}, 2u);
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(5.0f),  0.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_21_27_30) {
  PolynomApprox poly({32.0f, 48.0f, 63.0f}, {21.0f, 27.0f, 30.0f}, 2u);
  EXPECT_TRUE(eq(poly.eval(32.0f),  21.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(48.0f),  27.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(63.0f),  30.0f, 0.1f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f}, {1.0f, 4.0f, 9.0f}, 2u);
  EXPECT_TRUE(poly.size() == 3u);
  EXPECT_TRUE(eq(poly.eval(-1.0f), 1.0f));
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_TRUE(eq(poly.eval(4.0f), 16.0f));
  EXPECT_TRUE(eq(poly.eval(5.0f), 25.0f));
  EXPECT_TRUE(eq(poly.eval(6.0f), 36.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_5) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f}, {6.0f, 9.0f, 14.0f}, 2u);
  EXPECT_TRUE(poly.size() == 3u);
  EXPECT_TRUE(eq(poly.eval(-1.0f), 6.0f, 0.1));
  EXPECT_TRUE(eq(poly.eval(0.0f),  5.0f, 0.1));
  EXPECT_TRUE(eq(poly.eval(4.0f), 21.0f, 0.1));
  EXPECT_TRUE(eq(poly.eval(5.0f), 30.0f, 0.1));
  EXPECT_TRUE(eq(poly.eval(6.0f), 41.0f, 0.1));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2off) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f, 4.0f}, {1.0f, 4.0f, 9.0f, 15.9f}, 2u);
  EXPECT_TRUE(eq(poly.eval(0.0f),  0.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(4.0f), 16.0f, 0.1f));
  EXPECT_TRUE(eq(poly.eval(5.0f), 25.0f, 0.3f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.00059f));
}

TEST(polynomApprox, x2veryoff) {
  PolynomApprox poly({1.0f, 2.0f, 3.0f, 4.0f}, {1.0f, 3.0f, 13.0f, 13.0f}, 2u);
  EXPECT_FALSE(eq(poly.eval(0.0f),  0.0f));
  EXPECT_FALSE(eq(poly.eval(4.0f), 16.0f));
  EXPECT_FALSE(eq(poly.eval(5.0f), 25.0f));
  EXPECT_TRUE(poly.getRrmsError() > 0.1f);
}

TEST(polynomApprox, 2x3_3x2_4x_5) {
  PolynomApprox poly({0.0f, 1.0f, 2.0f, 3.0f}, {5.0f, 14.0f, 41.0f, 98.0f}, 3u);
  EXPECT_TRUE(poly.size() == 4u);
  EXPECT_TRUE(eq(poly.eval(0.5f),  8.0f));
  EXPECT_TRUE(eq(poly.eval(1.5f), 24.5f));
  EXPECT_TRUE(eq(poly.eval(2.5f), 65.0f));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(angle2apparentMirrorDepth, height2temp) {
  Angle2apparentMirrorDepth mirage(28);
  EXPECT_TRUE(eq(1.0, mirage.getHeightAtTemp(mirage.getTempAtHeight(1.0))));
  EXPECT_TRUE(eq(0.5, mirage.getHeightAtTemp(mirage.getTempAtHeight(0.5))));
  EXPECT_TRUE(eq(0.1, mirage.getHeightAtTemp(mirage.getTempAtHeight(0.1))));
  EXPECT_TRUE(eq(0.05, mirage.getHeightAtTemp(mirage.getTempAtHeight(0.05))));
  EXPECT_TRUE(eq(0.01, mirage.getHeightAtTemp(mirage.getTempAtHeight(0.01))));
  EXPECT_TRUE(eq(0.005, mirage.getHeightAtTemp(mirage.getTempAtHeight(0.005))));
}

TEST(angle2apparentMirrorDepth, temp2height) {
  Angle2apparentMirrorDepth mirage(28);
  EXPECT_TRUE(eq(300.0, mirage.getTempAtHeight(mirage.getHeightAtTemp(300.0))));
  EXPECT_TRUE(eq(308.0, mirage.getTempAtHeight(mirage.getHeightAtTemp(308.0))));
  EXPECT_TRUE(eq(316.0, mirage.getTempAtHeight(mirage.getHeightAtTemp(316.0))));
  EXPECT_TRUE(eq(324.0, mirage.getTempAtHeight(mirage.getHeightAtTemp(324.0))));
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
