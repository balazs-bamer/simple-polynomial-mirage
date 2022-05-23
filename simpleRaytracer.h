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
  Object(char const * const aName, Vertex const& aUpperCorner, Vertex const& aLowerCorner);
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
public:
  struct Parameters {
    uint32_t mRestrictCpu;
    double   mCamCenter;
    double   mTilt;
    double   mPinholeDist;
    double   mResolution;
    uint32_t mSubsample;
    uint32_t mGridColor;
    double   mGridIndent;
    uint32_t mGridSpacing;
  };

private:
  static constexpr double   csFilmSize        =  0.1;  // meters
  static constexpr double   csLimitHigh       =  cgPi / 33.3;
  static constexpr double   csLimitLow        = -cgPi / 33.3;
  static constexpr double   csLimitDelta      =  cgPi / 33333.3;
  static constexpr double   csLimitEpsilon    =  cgPi / 1234567.8;
  static constexpr double   csLimitAngleBoost =  1.001;
  static constexpr double   csSurfaceDistance =   1000; // meters
  static constexpr double   csSurfPinholeDist =      1; // meters
  static constexpr uint32_t csSurfSubsample   =      5u;

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
  double   const  mGridColor;
  double   const  mGridIndent;
  uint32_t const  mGridSpacing;
  Medium         &mMedium;
  double          mLimitAngleTop;
  double          mLimitAngleBottom;
  double          mLimitAngleDeep;
  double          mLimitAngleShallow;
  int             mLimitPixelBottom;
  int             mLimitPixelDeep;
  int             mLimitPixelShallow;

public:
  Image(Parameters const& aPara, Medium &aMedium);

  void process(char const * const aNameSurf, char const * const aNameOut);

private:
  void calculateLimits();
  void calculateMirage();
  void renderSurface(char const * const aNameSurf);

  static Vector getDirectionInXy(double const aAngle) { return Vector(std::cos(aAngle), std::sin(aAngle), 0.0); }
  static Vector getDirectionYz(double const aAngleY, double const aAngleZ) { return Vector(std::cos(aAngleY) * std::cos(aAngleZ), std::sin(aAngleY), std::cos(aAngleY) * std::sin(aAngleZ)); }
};

#endif
