#include "RungeKuttaRayBending.h"
#include "ShepardInterpolation.h"
#include "gtest/gtest.h"
#include <random>


constexpr float cgEpsilon = 0.0001f;

bool eq(double const aF1, double const aF2, double const aEpsilon) {
  return ::abs(aF1 - aF2) < aEpsilon;
}

bool eq(double const aF1, double const aF2) {
  return eq(aF1, aF2, cgEpsilon);
}

TEST(polynomApprox, x20) {
  double const y[] = {0.0, 0.0, 0.0};
  double const x[] = {1.0, 2.0, 3.0};
  PolynomApprox poly(3u, y, {{x, 2u}});
  EXPECT_TRUE(eq(poly.eval({0.0}),  0.0));
  EXPECT_TRUE(eq(poly.eval(4.0),  0.0));
  EXPECT_TRUE(eq(poly.eval(5.0),  0.0));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_21_27_30) {
  PolynomApprox poly(std::vector({21.0, 27.0, 30.0}), {{std::vector({32.0, 48.0, 63.0}), 2u}});
  EXPECT_TRUE(eq(poly.eval(32.0),  21.0, 0.1));
  EXPECT_TRUE(eq(poly.eval({48.0}),  27.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(63.0),  30.0, 0.1));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2) {
  PolynomApprox poly(std::vector({1.0, 4.0, 9.0}), {{std::vector({1.0, 2.0, 3.0}), 2u}});
  EXPECT_TRUE(eq(poly.eval(-1.0), 1.0));
  EXPECT_TRUE(eq(poly.eval(0.0),  0.0));
  EXPECT_TRUE(eq(poly.eval(4.0), 16.0));
  EXPECT_TRUE(eq(poly.eval(5.0), 25.0));
  EXPECT_TRUE(eq(poly.eval(6.0), 36.0));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2_5) {
  PolynomApprox poly(std::vector({6.0, 9.0, 14.0}), {{std::vector({1.0, 2.0, 3.0}), 2u}});
  EXPECT_TRUE(eq(poly.eval(-1.0), 6.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(0.0),  5.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(4.0), 21.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(5.0), 30.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(6.0), 41.0, 0.1));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, x2off) {
  PolynomApprox poly(std::vector({1.0, 4.0, 9.0, 15.9}), {{std::vector({1.0, 2.0, 3.0, 4.0}), 2u}});
  EXPECT_TRUE(eq(poly.eval(0.0),  0.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(4.0), 16.0, 0.1));
  EXPECT_TRUE(eq(poly.eval(5.0), 25.0, 0.3));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.00059f));
}

TEST(polynomApprox, x2veryoff) {
  PolynomApprox poly(std::vector({1.0, 3.0, 13.0, 13.0}), {{std::vector({1.0, 2.0, 3.0, 4.0}), 2u}});
  EXPECT_FALSE(eq(poly.eval(0.0),  0.0));
  EXPECT_FALSE(eq(poly.eval(4.0), 16.0));
  EXPECT_FALSE(eq(poly.eval(5.0), 25.0));
  EXPECT_TRUE(poly.getRrmsError() > 0.1f);
}

TEST(polynomApprox, 2x3_3x2_4x_5) {
  PolynomApprox poly(std::vector({5.0, 14.0, 41.0, 98.0}), {{std::vector({0.0, 1.0, 2.0, 3.0}), 3u}});
  EXPECT_TRUE(eq(poly.eval(0.5),  8.0));
  EXPECT_TRUE(eq(poly.eval(1.5), 24.5));
  EXPECT_TRUE(eq(poly.eval(2.5), 65.0));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, 10x_y_y2) {
  double const y[] = {0.0, 2.0, 6.0, 10.0, 12.0, 16.0, 20.0, 22.0, 26.0};
  double const x1[] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
  double const x2[] = {0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0};
  PolynomApprox poly(9u, y, {{x1, 1u}, {x2, 2u}});
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{0.0, 3.0}),  12.0));
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{3.0, 0.0}),  30.0));
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{3.0, 3.0}),  42.0));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, 10x_y_y2_xy) {
  double const y[] = {0.0, 2.0, 6.0, 10.0, 13.0, 18.0, 20.0, 24.0, 30.0};
  double const x1[] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
  double const x2[] = {0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0};
  PolynomApprox poly(9u, y, {{x1, 1u}, {x2, 2u}});
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{0.0, 3.0}),  12.0));
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{3.0, 0.0}),  30.0));
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{3.0, 3.0}),  51.0));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, 10y_x_x2_xy) {
  double const y[] = {0.0, 2.0, 6.0, 10.0, 13.0, 18.0, 20.0, 24.0, 30.0};
  double const x1[] = {0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0};
  double const x2[] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0};
  PolynomApprox poly(9u, y, {{x1, 2u}, {x2, 1u}});
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{0.0, 3.0}),  30.0));
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{3.0, 0.0}),  12.0));
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{3.0, 3.0}),  51.0));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

TEST(polynomApprox, 10y_x_x2__x3_xy) {
  double const y[] = {0.0, 1.0, -2.0, -15.0, 10.0, 12.0, 10.0, -2.0, 20.0, 23.0, 22.0, 11.0};
  double const x1[] = {0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0};
  double const x2[] = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0};
  PolynomApprox poly(9u, y, {{x1, 3u}, {x2, 1u}});
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{0.0, 4.0}),  40.0));
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{4.0, 0.0}), -44.0));
  EXPECT_TRUE(eq(poly.eval(std::initializer_list<double>{4.0, 4.0}),  12.0));
  EXPECT_TRUE(eq(poly.getRrmsError(), 0.0f));
}

/*TEST(shepardInterpolation, level0_dim1_data0) {
  using ShepIntpol = ShepardInterpolation<float, 1u, int, 4>;
  std::vector<ShepIntpol::Data> data;
  ShepIntpol shep(data, 3u);
  EXPECT_TRUE(shep.getTargetLevel() == std::numeric_limits<uint32_t>::max());
  EXPECT_TRUE(shep.getLevelCount() == 1u);
  EXPECT_TRUE(shep.getNodeCount(0u) == 1u);
  EXPECT_TRUE(shep.getItemCount(0u) == 0u);
}*/

TEST(fixedStack, iterator) {
  FixedStack<int, 10, 10> s;
  EXPECT_TRUE(s.isValid());
  for(int i = 0; i < 99; ++i) {
    s.push_back(i);
    EXPECT_TRUE(s.isValid());
  }
  EXPECT_TRUE(s.size() == 99u);
  for(int i = 98; i >= 0; --i) {
    EXPECT_TRUE(s.back() == i);
    s.pop_back();
    EXPECT_TRUE(s.isValid());
  }
  EXPECT_TRUE(s.size() == 0u);
}

TEST(fixedStack, clear) {
  FixedStack<int, 10, 10> s;
  for(int i = 0; i < 99; ++i) {
    s.push_back(i);
  }
  EXPECT_TRUE(s.size() == 99u);
  s.clear();
  EXPECT_TRUE(s.isValid());
  EXPECT_TRUE(s.size() == 0u);
}

TEST(shepardInterpolation, level10_dim1_data10) {
  using ShepIntpol = ShepardInterpolation<float, 1u, uint32_t, 3, 1>;
  typename ShepIntpol::DataTransfer data;
  for(uint32_t i = 0u; i < 10u; ++i) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = ::pow(1.0 / 3.0, i);
    item.mPayload = i;
    data.push_back(item);
  }
  ShepIntpol shep(data, 3u);
std::cout << "nodes: ";
for(uint32_t i = 0; i < shep.getLevelCount(); ++i) {
  std::cout << shep.getNodeCount(i) << ", ";
}
std::cout << "\nitems: ";
for(uint32_t i = 0; i < shep.getLevelCount(); ++i) {
  std::cout << shep.getItemCount(i) << ", ";
}
std::cout << "\n TL:" << shep.getTargetLevel() << '\n';
  EXPECT_TRUE(shep.getTargetLevel() == 8u);
  EXPECT_TRUE(shep.getLevelCount() == 10u);
  EXPECT_TRUE(shep.getNodeCount(0u) == 1u);
  EXPECT_TRUE(shep.getItemCount(0u) == 0u);
  EXPECT_TRUE(shep.getNodeCount(1u) == 2u);
  EXPECT_TRUE(shep.getItemCount(1u) == 1u);
  EXPECT_TRUE(shep.getNodeCount(2u) == 2u);
  EXPECT_TRUE(shep.getItemCount(2u) == 1u);
  EXPECT_TRUE(shep.getNodeCount(3u) == 1u);
  EXPECT_TRUE(shep.getItemCount(3u) == 0u);
  EXPECT_TRUE(shep.getNodeCount(4u) == 2u);
  EXPECT_TRUE(shep.getItemCount(4u) == 1u);
  EXPECT_TRUE(shep.getNodeCount(5u) == 2u);
  EXPECT_TRUE(shep.getItemCount(5u) == 1u);
  EXPECT_TRUE(shep.getNodeCount(6u) == 2u);
  EXPECT_TRUE(shep.getItemCount(6u) == 1u);
  EXPECT_TRUE(shep.getNodeCount(7u) == 2u);
  EXPECT_TRUE(shep.getItemCount(7u) == 1u);
  EXPECT_TRUE(shep.getNodeCount(8u) == 1u);
  EXPECT_TRUE(shep.getItemCount(8u) == 0u);
  EXPECT_TRUE(shep.getNodeCount(9u) == 2u);
  EXPECT_TRUE(shep.getItemCount(9u) == 4u);
  EXPECT_TRUE(shep.getDistanceFromTargetCenter(0.0) < 1.0);
  EXPECT_TRUE(shep.getDistanceFromTargetCenter(1.0) < 1.0);
}

TEST(shepardInterpolation, levelFew_dim1_random40) {
  using ShepIntpol = ShepardInterpolation<float, 1u, uint32_t, 3, 1>;
  typename ShepIntpol::DataTransfer data;
  uint32_t number = 0u;
  for(uint32_t i = 0u; i < 40u; ++i) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = number;
    item.mPayload = i;
    data.push_back(item);
    number = (number + 41u) % 64u;
  }
  ShepIntpol shep(data, 3u);
  EXPECT_TRUE(shep.getTargetLevel() == 3u);
  for(uint32_t i = 0u; i < 40u; ++i) {
    typename ShepIntpol::Location loc;
    loc[0] = number;
    EXPECT_TRUE(shep.getDistanceFromTargetCenter(loc) < 1.0);
    number = (number + 41u) % 64u;
  }
}

TEST(shepardInterpolation, levelFew_dim2_random40_data1d) {
  using ShepIntpol = ShepardInterpolation<float, 2u, uint32_t, 3, 1>;
  typename ShepIntpol::DataTransfer data;
  uint32_t number = 0u;
  for(uint32_t i = 0u; i < 40u; ++i) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = number;
    item.mLocation[1] = 0;
    item.mPayload = i;
    data.push_back(item);
    number = (number + 41u) % 64u;
  }
  ShepIntpol shep(data, 3u);
  EXPECT_TRUE(shep.getTargetLevel() == 3u);
}

TEST(shepardInterpolation, levelFew_dim2_random40_data2d) {
  using ShepIntpol = ShepardInterpolation<float, 2u, uint32_t, 3, 1>;
  typename ShepIntpol::DataTransfer data;
  uint32_t number = 0u;
  for(uint32_t i = 0u; i < 40u; ++i) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = number;
    number = (number + 43u) % 64u;
    item.mLocation[1] = number;
    number = (number + 43u) % 64u;
    item.mPayload = i;
    data.push_back(item);
  }
  ShepIntpol shep(data, 3u);
  EXPECT_TRUE(shep.getTargetLevel() == 2u);
  for(uint32_t i = 0u; i < 40u; ++i) {
    typename ShepIntpol::Location loc;
    loc[0] = number;
    number = (number + 41u) % 64u;
    loc[1] = number;
    number = (number + 41u) % 64u;
    EXPECT_TRUE(shep.getDistanceFromTargetCenter(loc) < 1.0);
  }
}

TEST(shepardInterpolation, levelFew_dim3_random40_data3d) {
  using ShepIntpol = ShepardInterpolation<float, 3u, CoefficientWise<double, 1u>, 3, 1 >;
  typename ShepIntpol::DataTransfer data;
  uint32_t number = 0u;
  for(uint32_t i = 0u; i < 97u; ++i) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = number;
    number = (number + 11u) % 41u;
    item.mLocation[1] = number;
    number = (number + 11u) % 41u;
    item.mLocation[2] = number;
    number = (number + 11u) % 41u;
    item.mPayload[0] = i;
    data.push_back(item);
  }
  ShepIntpol shep(data, 3u);
  EXPECT_TRUE(shep.getTargetLevel() == 3u);
  for(uint32_t i = 0u; i < 40u; ++i) {
    typename ShepIntpol::Location loc;
    loc[0] = number;
    number = (number + 13u) % 41u;
    loc[1] = number;
    number = (number + 13u) % 41u;
    loc[2] = number;
    number = (number + 13u) % 41u;
    shep.interpolate(loc);
    EXPECT_TRUE(shep.getDistanceFromTargetCenter(loc) < 1.0);
  }
}

/*TEST(angle2apparentMirrorDepth, height2temp) {
  ShepardRayBending mirage(28);
  EXPECT_TRUE(eq(1.0, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(1.0))));
  EXPECT_TRUE(eq(0.5, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.5))));
  EXPECT_TRUE(eq(0.1, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.1))));
  EXPECT_TRUE(eq(0.05, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.05))));
  EXPECT_TRUE(eq(0.01, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.01))));
  EXPECT_TRUE(eq(0.005, mirage.getHeightAtTempRise(mirage.getTempRiseAtHeight(0.005))));
}

TEST(angle2apparentMirrorDepth, temp2height) {
  ShepardRayBending mirage(28);
  EXPECT_TRUE(eq(0.0, mirage.getTempRiseAtHeight(mirage.getHeightAtTempRise(0.0))));
  EXPECT_TRUE(eq(8.0, mirage.getTempRiseAtHeight(mirage.getHeightAtTempRise(8.0))));
  EXPECT_TRUE(eq(16.0, mirage.getTempRiseAtHeight(mirage.getHeightAtTempRise(16.0))));
  EXPECT_TRUE(eq(24.0, mirage.getTempRiseAtHeight(mirage.getHeightAtTempRise(24.0))));
}*/

/* TODO test

  if(argc < 6) { return 1; }
  double bias = std::stod(std::string(argv[1]));
  std::cout << "bias: " << bias << '\n';
  double delta = std::stod(std::string(argv[2]));
  std::cout << "delta: " << delta << '\n';
  uint32_t toConsider = std::stoul(std::string(argv[3]));
  std::cout << "toConsider: " << toConsider << '\n';
  double shepardExponent = std::stod(std::string(argv[4]));
  std::cout << "shepExp: " << shepardExponent << '\n';
  double avRelSize = std::stod(std::string(argv[5]));
  std::cout << "avRelSize: " << avRelSize << '\n';
  using Data = CoefficientWise<double, 1u>;
  using ShepIntpol = ShepardInterpolation<double, 1u, Data, 4, 3>;
  std::vector<ShepIntpol::Data> data;
  double number = 0.0;
  double const max = 100.0;
  for(double number = 0.0; number < max; number += delta) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = number;
    // auto sum = number * number;          // relerr<0.007 if bias: 4321 delta: 2.3 toConsider: 4 shepExp: 3 TL: 3
    auto sum = std::sin(number / 7);        // biaserr<0.05 if bias: 2 delta: 2.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.006 unbiaserr<0.06 if bias: 12 delta: 2.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.003 unbiaserr<0.033 if bias: 22 delta: 2.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.002 unbiaserr<0.04 if bias: 42 delta: 2.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.002 unbiaserr<0.033 if bias: 22 delta: 1.3 toConsider: 4 shepExp: 3
                                            // biaserr<0.002 unbiaserr<0.033 if bias: 22 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 2
                                            // biaserr<0.0007 unbiaser<0.007 if bias: 22 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 3
                                            // biaserr<0.0007 unbiaser<0.007 if bias: 22 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 4
                                            // optimal bias = |vertspan| * 22
// biaserr is what seen on biased relative error graph
    item.mPayload = {static_cast<double>(sum)};
    data.push_back(item);
  }
  ShepIntpol shep(data, toConsider, avRelSize, shepardExponent, bias);
  std::cout << "TL: " << shep.getTargetLevel() << '\n';
  double diffSum = 0.0;
  double valSum = 0.0;
  double d2 = 0.25;
  std::cout << "diff=[";
  for(double i = 0u; i < max; i+=d2) {
    typename ShepIntpol::Location loc{i};
    // auto v = i*i;
    auto v = std::sin(i / 7);
    valSum += v;
    auto d = (v == 0.0 ? 1.0 : shep.interpolate(loc)[0] / v);
    diffSum += d * d;
    std::cout << d << (i < max - d2 ? ", " : "];\n");
  }
  auto diffAvg = std::sqrt(diffSum) / max / max;
  auto valAvg = valSum / max / max;
//  std::cout << "err: " << diffAvg/valAvg << '\n';


  using Data = CoefficientWise<double, 1u>;
  using ShepIntpol = ShepardInterpolation<double, 2u, Data, 6, 3>;
                                            // biaserr<0.0007 unbiaser<0.015 if bias: 17000 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 3
                                            // biaserr<0.0008 unbiaser<0.015 if bias: 15000 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 3
                                            // biaserr<0.0006 unbiaser<0.015 if bias: 19000 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 3
                                            // biaserr<0.0006 unbiaser<0.015 if bias: 19000 delta: 2.3 toConsider: 4 shepExp: 3 avgRelS: 0.25 avgCnt1d: 2
                                            // optimal bias = |vertspan| * 22
  std::vector<ShepIntpol::Data> data;
  double const max = 64.0;
  for(double n1 = 0.0; n1 < max; n1 += delta) {
    for(double n2 = 0.0; n2 < max; n2 += delta) {
      typename ShepIntpol::Data item;
      item.mLocation[0] = n1;
      item.mLocation[1] = n2 * 10.0;
      auto sum = (n1 * n1 + n2 * n2) / 10.0;
      item.mPayload = {sum};
      data.push_back(item);
    }
  }
  ShepIntpol shep(data, toConsider, avRelSize, shepardExponent, bias);
  std::cout << "TL: " << shep.getTargetLevel() << '\n';
  double diffSum = 0.0;
  double valSum = 0.0;
  std::cout << "diff=[";
  uint32_t limit = 64u;
  for(uint32_t i = 0u; i < limit; ++i) {
    for(uint32_t j = 0u; j < limit; ++j) {
      typename ShepIntpol::Location loc{i, j * 10.0};
      auto v = (i*i+j*j)/10.0;
      valSum += v;
      auto d = (v == 0.0 ? 1.0 : shep.interpolate(loc)[0] / v);
      diffSum += d * d;
      std::cout << d << (j < limit - 1u ? ", " : (i < limit - 1u ? ";\n" : "];\n"));
    }
  }
  auto diffAvg = std::sqrt(diffSum) / limit / limit;
  auto valAvg = valSum / limit / limit;
  std::cout << "err: " << diffAvg/valAvg << '\n';

  double slice = std::stod(std::string(argv[6]));
  std::cout << "slice: " << avRelSize << '\n';
  using Data = CoefficientWise<double, 1u>;
  using ShepIntpol = ShepardInterpolation<double, 3u, Data, 6, 3>;

  typename ShepIntpol::DataTransfer data;
  double const max = 64.0;
  for(double n1 = 0.0; n1 < max; n1 += delta) {
    for(double n2 = 0.0; n2 < max; n2 += delta) {
      for(double n3 = 0.0; n3 < max; n3 += delta) {
        typename ShepIntpol::Data item;
        item.mLocation[0] = n1;
        item.mLocation[1] = n2 *  5.0;
        item.mLocation[2] = n3 * 20.0;
        auto sum = (n1 * n1 + n2 * n2 + n3 * n3) / 10.0;
        item.mPayload = {sum};
        data.push_back(item);
      }
    }
  }
  ShepIntpol shep(data, toConsider, avRelSize, shepardExponent, bias);
  std::cout << "TL: " << shep.getTargetLevel() << '\n';
  double diffSum = 0.0;
  double valSum = 0.0;
  std::cout << "diff=[";
  uint32_t limit = 64u;
  for(uint32_t i = 0u; i < limit; ++i) {
    for(uint32_t j = 0u; j < limit; ++j) {
      typename ShepIntpol::Location loc{slice, i * 5.0, j * 20.0};
      auto v = (slice*slice+i*i+j*j)/10.0;
      valSum += v;
      auto d = (v == 0.0 ? 1.0 : shep.interpolate(loc)[0] / v);
      diffSum += d * d;
      std::cout << d << (j < limit - 1u ? ", " : (i < limit - 1u ? ";\n" : "];\n"));
    }
  }
  auto diffAvg = std::sqrt(diffSum) / limit / limit;
  auto valAvg = valSum / limit / limit;
  std::cout << "err: " << diffAvg/valAvg << '\n';


  using Data = CoefficientWise<double, 1u>;
  using ShepIntpol = ShepardInterpolation<double, 2u, Data, 3>;
  std::vector<ShepIntpol::Data> data;
  uint32_t number = 0u;
  for(uint32_t i = 0u; i < 60u; ++i) {
    typename ShepIntpol::Data item;
    item.mLocation[0] = number;
    auto sum = number * number;
    number = (number + 81u) % 127u;
    item.mLocation[1] = number;
    sum += number * number;
    number = (number + 81u) % 127u;
    item.mPayload = {static_cast<double>(sum)};
    data.push_back(item);
  }
  ShepIntpol shep(data, 3u, argc);
  std::cout << "TL: " << shep.getTargetLevel() << '\n';
  double diffSum = 0.0;
  double valSum = 0.0;
  for(uint32_t i = 0u; i < 127u; ++i) {
    for(uint32_t j = 0u; j < 127u; ++j) {
      typename ShepIntpol::Location loc{i, j};
      auto v = i*i+j*j;
      valSum += v;
      auto d = v - shep.interpolate(loc)[0];
      diffSum += d * d;
//    std::cout << "dftc: " << shep.getDistanceFromTargetCenter(loc) << " diffInt: " << d << '\n';
    }
  }
  auto diffAvg = ::sqrt(diffSum) / 127.0 / 127.0;
  auto valAvg = valSum / 127.0 / 127.0;
  std::cout << "err: " << diffAvg/valAvg << '\n';


  ShepardRayBending::Gather g;
  g.mAsphalt = false;
  g.mCollection.emplace_back(1.0, 2.0, 2.0/::sqrt(0.1*0.1+2.0*2.0));
  g.mCollection.emplace_back(0.9, 3.0, 3.0/::sqrt(0.1*0.1+3.0*3.0));
  g.mCollection.emplace_back(0.8, 4.0, 4.0/::sqrt(0.1*0.1+4.0*4.0));
  g.mCollection.emplace_back(0.75, 0.0, std::nan("1"));

  ShepardRayBending::Gather g;
  g.mAsphalt = false;
  g.mCollection.emplace_back(1.0, 2.0, 2.0/::sqrt(0.1*0.1+2.0*2.0));
  g.mCollection.emplace_back(0.9, 0.0, std::nan("1"));

  ShepardRayBending::Gather g;
  g.mAsphalt = true;
  g.mCollection.emplace_back(1.0, 2.0, 2.0/::sqrt(0.1*0.1+2.0*2.0));
  g.mCollection.emplace_back(0.9, 3.0, 3.0/::sqrt(0.1*0.1+3.0*3.0));
  g.mCollection.emplace_back(0.8, 4.0, 4.0/::sqrt(0.1*0.1+4.0*4.0));

  auto p = ShepardRayBending::toRayPath(g);
  std::cout << "d = [";
  for(auto i : p) {
    std::cout << i.mHorizDisp << ", ";
  }
  std::cout << "\nh = [";
  for(auto i : p) {
    std::cout << i.mHeight << ", ";
  }
  std::cout << "\na = [";
  for(auto i : p) {
    std::cout << i.mAngleFromHoriz << ", ";
  }
  std::cout << '\n';

  ShepardRayBending prb(33);
  auto p = prb.getRandomIndices(101,15);
  std::cout << "d = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(23,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(3,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(16,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(17,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(55,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(56,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(57,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(58,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(222,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  p = prb.getRandomIndices(323,15);
  std::cout << "\nd = [";
  for(auto i : p) {
    std::cout << i << ", ";
  }
  std::cout << '\n';
 */

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
