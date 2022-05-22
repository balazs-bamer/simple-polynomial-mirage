#include "simpleRaytracer.h"
#include "CLI11.hpp"
#include <iostream>
#include <thread>


int main(int aArgc, char **aArgv) {
  RungeKuttaRayBending::Parameters paraRk;
  Image::Parameters                paraIm;

  CLI::App opt{"Usage"};
  std::string nameBase = "water";
  opt.add_option("--base", nameBase, "base type (conventional / porous / water) [water]");
  double bullLift = 0.0;
  opt.add_option("--bullLift", bullLift, "lift of bulletin from ground (m) [0.0]");
  paraIm.mCamCenter = 1.1;
  opt.add_option("--camCenter", paraIm.mCamCenter, "height of camera center (m) [1.1]");
  double dist = 1000.0;
  opt.add_option("--dist", dist, "distance of bulletin and camera [1000]");
  std::string nameForm = "round";
  opt.add_option("--earthForm", nameForm, "Earth form (flat / round) [round]");
  double rawRadius = 6371.0;
  opt.add_option("--earthRadius", rawRadius, "Earth radius (km) [6371.0]");
  paraIm.mGridColor = 42u;
  opt.add_option("--gridColor", paraIm.mGridColor, "grid color, 0-black, 255-white (natural integer) [42]");
  paraIm.mGridIndent = 0.5;
  opt.add_option("--gridIndent", paraIm.mGridIndent, "grid indent ratio (0-1) from vertical edges compared to the black column (double) [0.5]");
  paraIm.mGridSpacing = 0u;
  opt.add_option("--gridSpacing", paraIm.mGridSpacing, "grid spacing in pixels, none if <2 (natural integer) [0]");
  double height = 9.0;
  opt.add_option("--height", height, "height of bulletin (m) [9.0]  its width will be calculated");
  paraRk.mMaxCosDirChange = 0.99999999999;
  opt.add_option("--maxCosDirChange", paraRk.mMaxCosDirChange, "Maximum of cos of direction change to reset big step [0.99999999999]");
  std::string nameIn = "monoscopeRca.png";
  opt.add_option("--nameIn", nameIn, "input filename [monoscopeRca.png]");
  std::string nameOut = "result.png";
  opt.add_option("--nameOut", nameOut, "output filename [result.png]");
  paraIm.mPinholeDist = 4.0;
  opt.add_option("--pinholeDist", paraIm.mPinholeDist, "pinhole distance from film (m) [4.0]");
  paraIm.mResolution = 1000u;
  opt.add_option("--resolution", paraIm.mResolution, "film resulution in both directions (pixel) [1000]");
  paraIm.mRestrictCpu = 0u;
  opt.add_option("--saveCpus", paraIm.mRestrictCpu, "amount of CPUs to save to keep the system responsive (natural integer) [0]");
  bool silent = false;
  opt.add_option("--silent", silent, "surpress parameter echo (true, false) [false]");
  paraRk.mStep1 = 0.01;
  opt.add_option("--step1", paraRk.mStep1, "initial step size (m) [0.01]");
  paraRk.mStepMin = 1e-4;
  opt.add_option("--stepMin", paraRk.mStepMin, "maximal step size (m) [1e-4]");
  paraRk.mStepMax = 55.5;
  opt.add_option("--stepMax", paraRk.mStepMax, "maximal step size (m) [55.5]");
  std::string nameStepper = "RungeKuttaFehlberg45";
  opt.add_option("--stepper", nameStepper, "stepper type (RungeKutta23 / RungeKuttaClass4 / RungeKuttaFehlberg45 / RungeKuttaCashKarp45 / RungeKuttaPrinceDormand89 / BulirschStoerBaderDeuflhard) [RungeKuttaFehlberg45]");
  paraIm.mSubsample = 2u;
  opt.add_option("--subsample", paraIm.mSubsample, "subsampling each pixel in both directions (count) [2]");
  double tempAmb = std::nan("");
  opt.add_option("--tempAmb", tempAmb, "ambient temperature (Celsius) [20 for conventional, 38.5 for porous, 10 for water]");
  double tempBase = 13.0;
  opt.add_option("--tempBase", tempBase, "base temperature, only for water (Celsius) [13]");
  paraIm.mTilt = 0.0;
  opt.add_option("--tilt", paraIm.mTilt, "camera tilt, neg downwards (degrees) [0.0]");
  paraRk.mTolAbs = 0.001;
  opt.add_option("--tolAbs", paraRk.mTolAbs, "absolute tolerance (m) [1e-3]");
  paraRk.mTolRel = 0.001;
  opt.add_option("--tolRel", paraRk.mTolRel, "relative tolerance (m) [1e-3]");
  CLI11_PARSE(opt, aArgc, aArgv);

  Eikonal::Model base;
  if(nameBase == "conventional") {
    base = Eikonal::Model::cConventional;
  }
  else if(nameBase == "porous") {
    base = Eikonal::Model::cPorous;
  }
  else if(nameBase == "water") {
    base = Eikonal::Model::cWater;
  }
  else {
    std::cerr << "Illegal base value: " << nameBase << '\n';
    return 1;
  }

  Eikonal::EarthForm earthForm;
  if(nameForm == "flat") {
    earthForm = Eikonal::EarthForm::cFlat;
  }
  else if(nameForm == "round") {
    earthForm = Eikonal::EarthForm::cRound;
  }
  else {
    std::cerr << "Illegal Earth form value: " << nameForm << '\n';
    return 1;
  }

  double earthRadius = rawRadius * 1000.0;

  if(nameStepper == "RungeKutta23") {
    paraRk.mStepper = StepperType::cRungeKutta23;
  }
  else if(nameStepper == "RungeKuttaClass4") {
    paraRk.mStepper = StepperType::cRungeKuttaClass4;
  }
  else if(nameStepper == "RungeKuttaFehlberg45") {
    paraRk.mStepper = StepperType::cRungeKuttaFehlberg45;
  }
  else if(nameStepper == "RungeKuttaCashKarp45") {
    paraRk.mStepper = StepperType::cRungeKuttaCashKarp45;
  }
  else if(nameStepper == "RungeKuttaPrinceDormand89") {
    paraRk.mStepper = StepperType::cRungeKuttaPrinceDormand89;
  }
  else if(nameStepper == "BulirschStoerBaderDeuflhard") {
    paraRk.mStepper = StepperType::cBulirschStoerBaderDeuflhard;
  }
  else {
    std::cerr << "Illegal stepper value: " << nameStepper << '\n';
    return 1;
  }

  if(std::isnan(tempAmb)) {
    tempAmb = (base == Eikonal::Model::cConventional ? 20.0 :
              (base == Eikonal::Model::cPorous ? 38.5 : 10.0));
  }
  else {} // nothing to do

  if(!silent) {
    std::cout << "base type:                                         " << nameBase << ' ' << static_cast<int>(base) << '\n';
    std::cout << "lift of bulletin from ground (m): .  .  .  .  .  . " << bullLift << '\n';
    std::cout << "height of camera center (m):                       " << paraIm.mCamCenter << '\n';
    std::cout << "distance of bulletin and camera (m):               " << dist << '\n';
    std::cout << "Earth form:        .  .  .  .  .  .  .  .  .  .  . " << nameForm << ' ' << static_cast<int>(earthForm) << '\n';
    std::cout << "Earth radius (km):                                 " << earthRadius / 1000.0 << '\n';
    std::cout << "grid color:                                        " << paraIm.mGridColor << '\n';
    std::cout << "grid indent: .  .  .  .  .  .  .  .  .  .  .  .  . " << paraIm.mGridIndent << '\n';
    std::cout << "grid spacing:                                      " << paraIm.mGridSpacing << '\n';
    std::cout << "height of bulletin (m):                            " << height << '\n';
    std::cout << "max of cos of direction change to reset big step:  " << std::setprecision(17) << paraRk.mMaxCosDirChange << '\n';
    std::cout << "input filename:                                    " << nameIn << '\n';
    std::cout << "output filename:                                   " << nameOut << '\n';
    std::cout << "pinhole distance from film (m):            .  .  . " << paraIm.mPinholeDist << '\n';
    std::cout << "film resolution in both directions (pixel):        " << paraIm.mResolution << '\n';
    std::cout << "initial step size (m):                             " << paraRk.mStep1 << '\n';
    std::cout << "minimal step size (m):   .  .  .  .  .  .  .  .  . " << paraRk.mStepMin << '\n';
    std::cout << "maximal step size (m):                             " << paraRk.mStepMax << '\n';
    std::cout << "stepper type:                                      " << nameStepper << ' ' << static_cast<int>(paraRk.mStepper) << '\n';
    std::cout << "subsampling each pixel in both directions (count): " << paraIm.mSubsample << '\n';
    std::cout << "ambient temperature (Celsius):                     " << tempAmb << '\n';
    std::cout << "base temperature, only for water (Celsius):        " << tempBase << '\n';
    std::cout << "camera tilt, neg downwards (degrees):      .  .  . " << paraIm.mTilt << '\n';
    std::cout << "absolute tolerance (m):                            " << paraRk.mTolAbs << '\n';
    std::cout << "relative tolerance (m):                            " << paraRk.mTolRel << '\n';
    uint32_t nCpus = std::thread::hardware_concurrency();
    nCpus -= (nCpus <= paraIm.mRestrictCpu ? nCpus - 1u : paraIm.mRestrictCpu);
    std::cout << "Using " << nCpus << " thread(s)" << std::endl;
  }
  else {} // nothing to do

  paraRk.mDistAlongRay    = dist * 2.0;
  auto effectiveRadius = (earthForm == Eikonal::EarthForm::cFlat ? std::numeric_limits<double>::infinity() : earthRadius);

  Object object(nameIn.c_str(), dist, bullLift, height, effectiveRadius);
  Medium medium(paraRk, earthForm, earthRadius, base, tempAmb, tempBase, object);
  Image image(paraIm, medium);
  image.calculateLimits();
  image.process(nameOut.c_str());
  return 0;
}
