#include  <boost/program_options.hpp>
#include "nr3.h"
#include "stepper.h"
#include "stepperdopr853.h"
#include "odeint.h"
#include "Eikonal.h"
#include "OdeSolver.h"
#include "OdeSolverGsl.h"
#include "StepperDormandPrice853.h"
#include "RungeKuttaRayBending.h"

//clang++ -Inumrec -I/usr/include/eigen3 -I../repos/eigen-initializer_list/src -DEIGEN_MATRIX_PLUGIN=\"Matrix_initializer_list.h\" -DEIGEN_ARRAY_PLUGIN=\"Array_initializer_list.h\" -std=c++17 eikonal.cpp mathUtil.cpp -o eikonal -ggdb -lboost_program_options


// These two are for conventional asphalt
double refract(double const aH, double const aTempAmb, double const aTempDiff) {
  auto celsius = aTempAmb + 0.018 + aTempDiff * std::exp(-aH * 10.08); // Toled vettem.
  return 1.0 + 7.86e-4 * 101 / (celsius + 273.15);
}

double refractDiff(double const aH, double const aTempAmb, double const aTempDiff) {
  auto t = aTempAmb + aTempDiff * std::exp(-10.08 * aH) + 273.168;
  return aTempDiff * 0.800211 * std::exp(-10.08 * aH) / t / t;
}


class Eik final {
public:
  static constexpr double cgC = 299792458; // m/2
  
private:
  double mTempAmb;
  double mTempDiff;

public:
  Eik(double const aTempAmb, double const aTempDiff) : mTempAmb(aTempAmb), mTempDiff(aTempDiff) {}

  void operator() (const double/* aS*/, VecDoub_I const &aY, VecDoub_O &aDyds) {
    double n    = refract(aY[2], mTempAmb, mTempDiff);
    double v    = cgC / n;
//    double dvdz = -v / n * refractDiff(aY[2]);

    aDyds[0] = v * aY[3];
    aDyds[1] = v * aY[4];
    aDyds[2] = v * aY[5];
    aDyds[3] = 0.0;
    aDyds[4] = 0.0;
    aDyds[5] = refractDiff(aY[2], mTempAmb, mTempDiff) / v / n;
    //aDyds[5] = -1.0 / v / v * dvdz;
  }
};

template<typename tReal>
class Eikonal2 final {
public:
  static constexpr double cgC = 299792458; // m/2
  static constexpr uint32_t csNvar = 6u;
  using Real                       = tReal;
  using Variables                  = std::array<Real, csNvar>;
  
private:
  double mTempAmb;
  double mTempDiff;

public:
  Eikonal2(double const aTempAmb, double const aTempDiff) : mTempAmb(aTempAmb), mTempDiff(aTempDiff) {}

  void operator() (const double/* aS*/, std::array<Real, csNvar> const &aY, std::array<Real, csNvar> &aDyds) const {
    double n    = refract(aY[2], mTempAmb, mTempDiff);
    double v    = cgC / n;
//    double dvdz = -v / n * refractDiff(aY[2]);

    aDyds[0] = v * aY[3];
    aDyds[1] = v * aY[4];
    aDyds[2] = v * aY[5];
    aDyds[3] = 0.0;
    aDyds[4] = 0.0;
    aDyds[5] = refractDiff(aY[2], mTempAmb, mTempDiff) / v / n;
    //aDyds[5] = -1.0 / v / v * dvdz;
  }
};

void comp(double aDir, double aDist, double aHeight, double aStep1, double aStepMin, double aTempAmb, double aTempDiff, double aTolAbs, double aTolRel, double aTarget) {
  constexpr uint32_t cNvar = 6u;
  VecDoub yStart(cNvar);
  yStart[0] = 0.0;
  yStart[1] = 0.0;
  yStart[2] = aHeight;
  auto u = refract(aHeight, aTempAmb, aTempDiff) / Eik::cgC;
  yStart[3] = u * std::cos(aDir / 180.0 * 3.1415926539);
  yStart[4] = 0.0;
  yStart[5] = u * std::sin(aDir / 180.0 * 3.1415926539);
  Output out(100);
  Eik eik(aTempAmb, aTempDiff);
  Odeint<StepperDopr853<Eik>> ode(yStart, 0.0, aDist, aTolAbs, aTolRel, aStep1, aStepMin, out, eik);
  ode.integrate();
  std::cout << "xnro=[";
  for (Int i=0; i<out.count; i++) {
    std::cout << out.ysave[0][i] << (i < out.count - 1 ? ", " : "];\n");
  }
  std::cout << "znro=[";
  for (Int i=0; i<out.count; i++) {
    std::cout << out.ysave[2][i] << (i < out.count - 1 ? ", " : "];\n");
  }
  std::cout << "\n";

  std::array<double, cNvar> yStart2;
  yStart2[0] = 0.0;
  yStart2[1] = 0.0;
  yStart2[2] = aHeight;
  yStart2[3] = u * std::cos(aDir / 180.0 * 3.1415926539);
  yStart2[4] = 0.0;
  yStart2[5] = u * std::sin(aDir / 180.0 * 3.1415926539);
  Eikonal2<double> eikonal2(aTempAmb, aTempDiff);
  OdeSolver<StepperDormandPrice853<Eikonal2<double>>> ode2(0.0, aDist, aTolAbs, aTolRel, aStep1, aStepMin, eikonal2);
  auto judge = [&aTarget](std::array<double, cNvar> const& aY){ return aY[0] >= aTarget; };
  auto result = ode2.solve(yStart2, judge, 0.0001);
  std::cout << "xnrm=[0, " << result[0] << "];\n";
  std::cout << "znrm=[" << aHeight << ", " << result[2] << "];\n";
  aTarget /= 2;
  result = ode2.solve(yStart2, judge, 0.0001);
  std::cout << "xnrm2=[0, " << result[0] << "];\n";
  std::cout << "znrm2=[" << aHeight << ", " << result[2] << "];\n";

  Eikonal eikonal(aTempAmb, aTempDiff);
  RungeKuttaRayBending rk(aDist, aTolAbs, aTolRel, aStep1, aStepMin, eikonal);
  Vertex start(0.0f, 0.0f, (float)(aHeight));
  Vector dir((float)(std::cos(aDir / 180.0 * 3.1415926539)), 0.0f, (float)(std::sin(aDir / 180.0 * 3.1415926539)));
  auto final = rk.solve4x(start, dir, aTarget * 0.75);
  std::cout << "xnrrb=[0, " << final(0) << "];\n";
  std::cout << "znrrb=[" << aHeight << ", " << final(2) << "];\n";

  OdeSolverGsl<Eikonal> gsl(0.0, aDist, 1e-6, 1e-6, 1e-6, eikonal);
  typename Eikonal::Variables yStart3;
  yStart3[0] = 0.0;
  yStart3[1] = aHeight;
  yStart3[2] = 0.0;
  yStart3[3] = u * std::cos(aDir / 180.0 * 3.1415926539);
  yStart3[4] = u * std::sin(aDir / 180.0 * 3.1415926539);
  yStart3[5] = 0.0;
  std::vector<typename Eikonal::Variables> stuff;
  for(double t = 0.0; t < aDist; t += aDist / 100.0) {
    auto [tt,y] = gsl.solve(yStart3, [t](typename Eikonal::Variables const &aY){ return aY[0] > t; });
    stuff.push_back(y);
  }
  std::cout << "xgls=[";
  for (int i = 0; i < stuff.size(); ++i) {
    std::cout << stuff[i][0] << (i < stuff.size() - 1 ? ", " : "];\n");
  }
  std::cout << "zgls=[";
  for (int i = 0; i < stuff.size(); ++i) {
    std::cout << stuff[i][1] << (i < stuff.size() - 1 ? ", " : "];\n");
  }
  std::cout << "\n";
}

int main(int aArgc, char **aArgv) {
  boost::program_options::options_description desc("Usage:");
  desc.add_options()
                    ("dir", boost::program_options::value<double>(), "start direction (degrees), neg is down, 0 horizontal [0]")
                    ("dist", boost::program_options::value<double>(), "horizontal distance to track [2000]")
                    ("height", boost::program_options::value<double>(), "start height (m) [1]")
                    ("help", "produce this message")
                    ("step1", boost::program_options::value<double>(), "initial step size [0.01]")
                    ("stepmin", boost::program_options::value<double>(), "minimal step size [0.0]")
                    ("tempamb", boost::program_options::value<double>(), "ambient temperature (Celsius) [20]")
                    ("tempdiff", boost::program_options::value<double>(), "temperature rise next to asphalt compared to ambient (Celsius) [6.37]")
                    ("tolabs", boost::program_options::value<double>(), "absolute tolerance [1e-3]")
                    ("tolrel", boost::program_options::value<double>(), "relative tolerance [1e-3]")
                    ("xyz", boost::program_options::value<double>(), "custom [0]");
  boost::program_options::variables_map varMap;
  boost::program_options::store(boost::program_options::parse_command_line(aArgc, aArgv, desc), varMap);
  boost::program_options::notify(varMap);
  if (varMap.count("help")) {
    cout << desc << "\n";
  }
  else {
    double dir = (varMap.count("dir") ? varMap["dir"].as<double>() : 0.0);
    double dist = (varMap.count("dist") ? varMap["dist"].as<double>() : 2000.0);
    double height = (varMap.count("height") ? varMap["height"].as<double>() : 1.0);
    double step1 = (varMap.count("step1") ? varMap["step1"].as<double>() : 0.01);
    double stepMin = (varMap.count("stepmin") ? varMap["stepmin"].as<double>() : 0.0);
    double tempAmb = (varMap.count("tempamb") ? varMap["tempamb"].as<double>() : 20);
    double tempDiff = (varMap.count("tempdiff") ? varMap["tempdiff"].as<double>() : 6.37);
    double tolAbs = (varMap.count("tolabs") ? varMap["tolabs"].as<double>() : 0.001);
    double tolRel = (varMap.count("tolrel") ? varMap["tolrel"].as<double>() : 0.001);
    double target = (varMap.count("xyz") ? varMap["xyz"].as<double>() : 0.0);
    comp(dir, dist, height, step1, stepMin, tempAmb, tempDiff, tolAbs, tolRel, target);
  }
  return 0;
}
