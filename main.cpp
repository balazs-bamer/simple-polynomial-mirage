#include "simpleRaytracer.h"
#include "CLI11.hpp"
#include <iostream>


int main(int aArgc, char **aArgv) {
  CLI::App opt{"Usage"};
  double bullCenter = 5.0;
  opt.add_option("--bullCenter", bullCenter, "height of bulletin center (m) [5.0]");
  double camCenter = 1.1;
  opt.add_option("--camCenter", camCenter, "height of camera center (m) [1.1]");
  double dist = 1000.0;
  opt.add_option("--dist", dist, "distance of bulletin and camera [1000]");
  double height = 9.0;
  opt.add_option("--height", height, "height of bulletin (m) [9.0]");
  std::string nameIn = "monoscope.png";
  opt.add_option("--nameIn", nameIn, "input filename [monoscope.png]");
  std::string nameOut = "result.png";
  opt.add_option("--nameOut", nameOut, "output filename [result.png]");
  double pinholeDist = 4.0;
  opt.add_option("--pinholeDist", pinholeDist, "pinhole distance from film (m) [20.0]");
  uint32_t resolution = 1000u;
  opt.add_option("--resolution", resolution, "film resulution in both directions (pixel) [1000]");
  double step1 = 0.01;
  opt.add_option("--step1", step1, "initial step size [0.01]");
  uint32_t subsample = 2u;
  opt.add_option("--subsample", subsample, "subsampling each pixel in both directions (count) [2]");
  double tempAmb = 20.0;
  opt.add_option("--tempAmb", tempAmb, "ambient temperature (Celsius) [20]");
  double tempDiff = 6.37;
  opt.add_option("--tempDiff", tempDiff, "temperature rise next to asphalt compared to ambient (Celsius) [6.37]");
  double tilt = 0.0;
  opt.add_option("--tilt", tilt, "camera tilt in degrees, neg downwards [0.0]");
  double tolAbs = 0.001;
  opt.add_option("--tolAbs", tolAbs, "absolute tolerance [1e-3]");
  double tolRel = 0.001;
  opt.add_option("--tolRel", tolRel, "relative tolerance [1e-3]");
  double width = 12.0;
  opt.add_option("--width", width, "width of bulletin [12.0]");
  CLI11_PARSE(opt, aArgc, aArgv);

  Object object(nameIn.c_str(), dist, bullCenter, width, height);
  Medium medium(tempAmb, tempDiff, dist * 2.0, tolAbs, tolRel, step1, object);
  Image image(camCenter, tilt, pinholeDist, 0.1 / resolution, resolution, resolution, subsample, medium);
  image.process(nameOut.c_str());
  return 0;
}
