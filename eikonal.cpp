#include  <boost/program_options.hpp>
#include <iostream>
#include <iomanip>
#include "nr3.h"
#include "stepper.h"
#include "stepperdopr853.h"
#include "odeint.h"
#include "Eikonal.h"
#include "OdeSolver.h"
#include "OdeSolverGsl.h"
#include "StepperDormandPrice853.h"
#include "RungeKuttaRayBending.h"

//clang++ -Inumrec -I/usr/include/eigen3 -I../repos/eigen-initializer_list/src -DEIGEN_MATRIX_PLUGIN=\"Matrix_initializer_list.h\" -DEIGEN_ARRAY_PLUGIN=\"Array_initializer_list.h\" -std=c++17 eikonal.cpp mathUtil.cpp -o eikonal -ggdb -lboost_program_options -lgsl


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

void comp(StepperType aStepper, double aDir, double aDist, double aHeight, double aStep1, double aStepMin, double aTempAmb, Eikonal::Mode aMode, double aTolAbs, double aTolRel, double aTarget, uint32_t aSamples) {
/*  constexpr uint32_t cNvar = 6u;
  VecDoub yStart(cNvar);
  yStart[0] = 0.0;
  yStart[1] = 0.0;
  yStart[2] = aHeight;
  auto u = refract(aHeight, aTempAmb, 6.37) / Eik::cgC;
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
  std::cout << "znrm2=[" << aHeight << ", " << result[2] << "];\n";*/

  Eikonal eikonal(aTempAmb, aMode);
/*  RungeKuttaRayBending rk(aDist, aTolAbs, aTolRel, aStep1, aStepMin, eikonal);
  Vertex start(0.0f, 0.0f, (float)(aHeight));
  Vector dir((float)(std::cos(aDir / 180.0 * 3.1415926539)), 0.0f, (float)(std::sin(aDir / 180.0 * 3.1415926539)));
  auto final = rk.solve4x(start, dir, aTarget * 0.75);
  std::cout << "xnrrb=[0, " << final(0) << "];\n";
  std::cout << "znrrb=[" << aHeight << ", " << final(2) << "];\n";*/

  auto u = eikonal.getSlowness(aHeight);
  OdeSolverGsl<Eikonal> gsl(aStepper, 0.0, aDist, aTolAbs, aTolRel, aStep1, eikonal);
  typename Eikonal::Variables yStart3;
  yStart3[0] = 0.0;
  yStart3[1] = aHeight;
  yStart3[2] = 0.0;
  yStart3[3] = u * std::cos(aDir / 180.0 * 3.1415926539);
  yStart3[4] = u * std::sin(aDir / 180.0 * 3.1415926539);
  yStart3[5] = 0.0;
  std::vector<typename Eikonal::Variables> stuff;
  ofstream out("values.txt");
  for(double t = 0.0; t < aTarget; t += aTarget / aSamples) {
    auto [tt,y] = gsl.solve(yStart3, [t](typename Eikonal::Variables const &aY){ return aY[0] > t; });
    stuff.push_back(y);
    out << std::setprecision(10) << y[0] << '\t' << std::setprecision(10) << y[1] << '\n';
  }
  std::cout << "xgls=[";
  for (int i = 0; i < stuff.size(); ++i) {
    std::cout << std::setprecision(10) << stuff[i][0] << (i < stuff.size() - 1 ? ", " : "];\n");
  }
  std::cout << "zgls=[";
  for (int i = 0; i < stuff.size(); ++i) {
    std::cout << std::setprecision(10) << stuff[i][1] << (i < stuff.size() - 1 ? ", " : "];\n");
  }
  std::cout << "\n";
/*  std::cout << "h=[";
  for (auto i = 0.01; i < 1; i += .01) {
    std::cout << std::setprecision(10) << i << ", ";
  }
  std::cout << "n=[";
  for (auto i = 0.01; i < 1; i += .01) {
    std::cout << std::setprecision(10) << eikonal.getRefract(i) << ", ";
  }
  std::cout << "\n";*/
}

int main(int aArgc, char **aArgv) {
  boost::program_options::options_description desc("Usage:");
  desc.add_options()
                    ("asphalt", boost::program_options::value<std::string>(), "asphalt type (conventional / porous) [conventional]")
                    ("dir", boost::program_options::value<double>(), "start direction (degrees), neg is down, 0 horizontal [0]")
                    ("dist", boost::program_options::value<double>(), "distance along the ray to track [2000]")
                    ("height", boost::program_options::value<double>(), "start height (m) [1]")
                    ("help", "produce this message")
                    ("samples", boost::program_options::value<uint32_t>(), "number of samples on ray [100]")
                    ("step1", boost::program_options::value<double>(), "initial step size [1e-6]")
                    ("stepper", boost::program_options::value<std::string>(), "stepper type (RungeKutta23 / RungeKuttaClass4 / RungeKuttaFehlberg45 / RungeKuttaCashKarp45 / RungeKuttaPrinceDormand89 / BulirschStoerBaderDeuflhard) [RungeKutta23]")
                    ("stepmin", boost::program_options::value<double>(), "minimal step size [0.0]")
                    ("target", boost::program_options::value<double>(), "horizontal distance to travel [1000]")
                    ("tempamb", boost::program_options::value<double>(), "ambient temperature (Celsius) [20 for conventional, 38.5 for porous]")
                    ("tolabs", boost::program_options::value<double>(), "absolute tolerance [1e-6]")
                    ("tolrel", boost::program_options::value<double>(), "relative tolerance [1e-6]");
  boost::program_options::variables_map varMap;
  boost::program_options::store(boost::program_options::parse_command_line(aArgc, aArgv, desc), varMap);
  boost::program_options::notify(varMap);
  if (varMap.count("help")) {
    cout << desc << "\n";
  }
  else {
    Eikonal::Mode asphalt = Eikonal::Mode::cConventional;
    std::string name = (varMap.count("asphalt") ? varMap["asphalt"].as<std::string>() : std::string("conventional"));
    if(name == "porous") {
      asphalt = Eikonal::Mode::cPorous;
    }
    else {} // nothing to do
    double dir = (varMap.count("dir") ? varMap["dir"].as<double>() : 0.0);
    double dist = (varMap.count("dist") ? varMap["dist"].as<double>() : 2000.0);
    double height = (varMap.count("height") ? varMap["height"].as<double>() : 1.0);
    uint32_t samples = (varMap.count("samples") ? varMap["samples"].as<uint32_t>() : 100);
    double step1 = (varMap.count("step1") ? varMap["step1"].as<double>() : 1e-6);
    name = (varMap.count("stepper") ? varMap["stepper"].as<std::string>() : std::string("RungeKutta23"));
    StepperType stepper = StepperType::cRungeKutta23;
    if(name == "RungeKutta23") {
      stepper = StepperType::cRungeKutta23;
    }
    else if(name == "RungeKuttaClass4") {
      stepper = StepperType::cRungeKuttaClass4;
    }
    else if(name == "RungeKuttaFehlberg45") {
      stepper = StepperType::cRungeKuttaFehlberg45;
    }
    else if(name == "RungeKuttaCashKarp45") {
      stepper = StepperType::cRungeKuttaCashKarp45;
    }
    else if(name == "RungeKuttaPrinceDormand89") {
      stepper = StepperType::cRungeKuttaPrinceDormand89;
    }
    else if(name == "BulirschStoerBaderDeuflhard") {
      stepper = StepperType::cBulirschStoerBaderDeuflhard;
    }
    else {} // nothing to do
    double stepMin = (varMap.count("stepmin") ? varMap["stepmin"].as<double>() : 0.0);
    double target = (varMap.count("target") ? varMap["target"].as<double>() : 1000.0);
    double tempAmb = (varMap.count("tempamb") ? varMap["tempamb"].as<double>() : (asphalt == Eikonal::Mode::cConventional ? 20.0 : 38.5));
    double tempDiff = (varMap.count("tempdiff") ? varMap["tempdiff"].as<double>() : 6.37);
    double tolAbs = (varMap.count("tolabs") ? varMap["tolabs"].as<double>() : 1e-6);
    double tolRel = (varMap.count("tolrel") ? varMap["tolrel"].as<double>() : 1e-6);
    comp(stepper, dir, dist, height, step1, stepMin, tempAmb, asphalt, tolAbs, tolRel, target, samples);
  }
  return 0;
}
