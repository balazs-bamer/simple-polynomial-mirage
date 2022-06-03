#ifndef SIMPLERAYTRACER_H
#define SIMPLERAYTRACER_H

#define __FreeBSD__ 12 // Hack to let png++ compile under cygwin

#include "RungeKuttaRayBending.h"
#include "3dGeomUtil.h"
#include "png.hpp"
#include <optional>


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
         double const aTempAmbient, double const tempAmbMin, double const tempAmbMax, double const aTempBase, Object const& aObject)
  : mEikonal(aEarthForm, aEarthRadius, aModel, aTempAmbient, tempAmbMin, tempAmbMax, aTempBase)
  , mSolver(aParameters, mEikonal)
  , mObject(aObject) {}

  Medium(Medium const&) = default;
  Medium(Medium &&) = delete;
  Medium& operator=(Medium const&) = delete;
  Medium& operator=(Medium &&) = delete;

  void setWaterTempAmb(Eikonal::Temperature const aWhich) { mEikonal.setWaterTempAmb(aWhich); }
  uint8_t trace(Ray const& aRay);
  bool hits(Ray const& aRay);
  RungeKuttaRayBending::Result getHit(Ray const& aRay) { return mSolver.solve4x(aRay.mStart, aRay.mDirection, mObject.getX()); }
  double getRefract(double const aH) const { return mSolver.getRefract(aH); }
};


class Image final {
public:
  struct Parameters {
    uint32_t mRestrictCpu;
    double   mCamCenter;
    double   mTilt;
    double   mBorderFactor;
    uint32_t mResolutionX;
    uint32_t mSubsample;
    double   mMarkIndent;
    bool     mMarkAcross;
  };

private:
  static constexpr double   csLimitHigh           =  cgPi / 33.3;
  static constexpr double   csLimitLow            = -cgPi / 33.3;
  static constexpr double   csLimitDelta          =  cgPi / 33333.3;
  static constexpr double   csLimitEpsilon        =  cgPi / 1234567.8;
  static constexpr double   csLimitAngleBoost     =      1.001;
  static constexpr double   csRenderSurfaceFactor =      2.0;
  static constexpr double   csSurfaceDistance     =   1000; // meters
  static constexpr double   csSurfPinholeDist     =      1; // meters
  static constexpr uint32_t csSurfSubsample       =      5u;
  static constexpr uint8_t  csColorVoid           =      0u;
  static constexpr uint8_t  csColorMirror         =      1u;
  static constexpr uint8_t  csColorBase           =      2u;
  static constexpr uint8_t  csColorBlack          =      3u;
  static constexpr int      csDashCount           =     20;

  uint32_t const  mRestrictCpu;
  std::vector<uint8_t>         mBuffer;
  png::image<png::index_pixel> mImage;
  png::palette                 mPalette;
  uint32_t const  mSubSample;
  uint32_t const  mResolutionX;
  double   const  mBorderFactor;
  double   const  mSsFactor;
  Vertex   const  mCenter;
  Vector   const  mNormal;
  Vector   const  mInPlaneZ;
  Vector   const  mInPlaneY;
  Vertex   const  mPinhole;
  double   const  mBiasSub;
  double   const  mMarkIndent;
  bool     const  mMarkAcross;

  Medium                &mMedium;
  std::optional<double>  mLimitAngleTop;
  std::optional<double>  mLimitAngleBottom;
  std::optional<double>  mLimitAngleDeep;
  std::optional<double>  mLimitAngleShallow;
  double                 mLimitAngleBottomSurf;
  int                    mLimitPixelBaseTop;
  int                    mLimitPixelBaseBottom;
  int                    mLimitPixelBaseBottomSurf;
  int                    mLimitPixelDeep;
  int                    mLimitPixelShallow;
  int                    mLimitPixelBottom;
  double                 mPixelSize;
  double                 mBiasZ;
  double                 mBiasY;

public:
  Image(Parameters const& aPara, Medium &aMedium);

  void process(char const * const aNameSurf, char const * const aNameOut);

private:
  void calculateAngleLimits(Eikonal::Temperature const aWhich);
  void calculateBiases(bool const aRenderSurface);
  int calculatePixelLimitY(double const aAngle);
  int calculatePixelLimitZ(double const aAngle);
  int calculateMirrorHeight();
  void renderSurface(char const * const aNameSurf);
  void calculateMirage();
  void drawMarks(int const aMirrorHeight);

  static Vector getDirectionInXy(double const aAngle) { return Vector(std::cos(aAngle), std::sin(aAngle), 0.0); }
  static Vector getDirectionInXz(double const aAngle) { return Vector(std::cos(aAngle), 0.0, std::sin(aAngle)); }
  static Vector getDirectionYz(double const aAngleY, double const aAngleZ) { return Vector(std::cos(aAngleY) * std::cos(aAngleZ), std::sin(aAngleY), std::cos(aAngleY) * std::sin(aAngleZ)); }
};

#endif
