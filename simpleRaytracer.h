#ifndef SIMPLERAYTRACER_H
#define SIMPLERAYTRACER_H

#include "RungeKuttaRayBending.h"
#include "3dGeomUtil.h"
#include "png.hpp"

class Object final {
private:
  png::image<png::gray_pixel> mImage;
  double const mDy;
  double const mDz;
  double       mMinY;
  double       mMaxY;
  double const mMinZ;
  double const mMaxZ;
  double const mX;

public:
  Object(char const * const aName, double const aDispX, double const aLiftY, double const aHeight, double const aEarthRadius);
  double  getX() const { return mX; }
  bool    hasPixel(Vertex const &aHit) const;
  uint8_t getPixel(Vertex const &aHit) const;
};


class Medium final {
private:
  Eikonal              mEikonal;
  RungeKuttaRayBending mSolver;
  Object const&        mObject;

public:
  Medium(RungeKuttaRayBending::Parameters const& aParameters,
         Eikonal::EarthForm const aEarthForm, double const aEarthRadius, Eikonal::Model const aModel,
         double const aTempAmbient, double const aTempBase, Object const& aObject)
  : mEikonal(aEarthForm, aEarthRadius, aModel, aTempAmbient, aTempBase)
  , mSolver(aParameters, mEikonal)
  , mObject(aObject) {}

  Medium(Medium const&) = default;
  Medium(Medium &&) = delete;
  Medium& operator=(Medium const&) = delete;
  Medium& operator=(Medium &&) = delete;

  uint8_t trace(Ray const& aRay);
  bool hits(Ray const& aRay);
};


class Image final {
private:
  static constexpr double csLimitHigh    =  cgPi / 33.3;
  static constexpr double csLimitLow     = -cgPi / 33.3;
  static constexpr double csLimitDelta   =  cgPi / 33333.3;
  static constexpr double csLimitEpsilon =  cgPi / 1234567.8;

  uint32_t const  mRestrictCpu;
  std::vector<uint8_t>        mBuffer;
  png::image<png::gray_pixel> mImage;
  uint32_t const  mSubSample;
  double   const  mSsFactor;
  double   const  mPixelSize;
  Vertex   const  mCenter;
  Vector   const  mNormal;
  Vector   const  mInPlaneZ;
  Vector   const  mInPlaneY;
  Vertex   const  mPinhole;
  double   const  mBiasZ;
  double   const  mBiasY;
  double   const  mBiasSub;
  Medium         &mMedium;

public:
  Image(uint32_t const aRestrictCpu, double const aCenterY,
        double const aTilt, double const aPinholeDist,
        double const aPixelSize,
        uint32_t const aResZ, uint32_t const aResY,
        uint32_t const aSubSample, Medium &aMedium);

  void dumpLimits() const;
  void process(char const * const aName);

private:
  static Vector getDirectionInXy(double const aAngle) { return Vector(std::cos(aAngle), std::sin(aAngle), 0.0); }
};

#endif
