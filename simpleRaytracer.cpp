#include "simpleRaytracer.h"
#include <iomanip>
#include <iostream>
#include <thread>


Object::Object(char const * const aName, double const aDispX, double const aLiftY, double const aHeight, double const aEarthRadius)
  : mImage(aName)
  , mDy(aHeight / mImage.get_height())
  , mDz(mDy)
  , mMinY(aLiftY)
  , mMaxY(aLiftY + aHeight)
  , mMinZ(-static_cast<double>(mImage.get_width()) * aHeight / static_cast<double>(mImage.get_height()) / 2.0)
  , mMaxZ(-mMinZ)
  , mX(aDispX) {
  double shift = (std::isinf(aEarthRadius) ? 0.0 : std::sqrt(aEarthRadius * aEarthRadius - mX * mX) - aEarthRadius);
  mMinY += shift;
  mMaxY += shift;
}

Object::Object(char const * const aName, Vertex const& aUpperCorner, Vertex const& aLowerCorner)
  : mImage(aName)
  , mDy((aUpperCorner(1) - aLowerCorner(1)) / mImage.get_height())
  , mDz((aUpperCorner(2) - aLowerCorner(2)) / mImage.get_width())
  , mMinY(aLowerCorner(1))
  , mMaxY(aUpperCorner(1))
  , mMinZ(aLowerCorner(2))
  , mMaxZ(aUpperCorner(2))
  , mX(aLowerCorner(0)) {
}

bool Object::hasPixel(Vertex const &aHit) const {
  return aHit(1) > mMinY && aHit(1)  < mMaxY && aHit(2) > mMinZ && aHit(2) < mMaxZ;
}

uint8_t Object::getPixel(Vertex const &aHit) const {
  uint8_t result = 0u;
  if(hasPixel(aHit)) {
    uint32_t x = std::max(0.0, (aHit(2) - mMinZ) / mDz);
    uint32_t y = mImage.get_height() - (aHit(1) - mMinY) / mDy - 1u;
    result = mImage.get_pixel(x, y);
  }
  else {} // Nothing to do
  return result;
}


uint8_t Medium::trace(Ray const& aRay) {
  try {
    auto hit = mSolver.solve4x(aRay.mStart, aRay.mDirection, mObject.getX());
    if(hit.mValid) {
      return mObject.getPixel(hit.mValue);
    }
    else {
      return 0;
    }
  }
  catch(...) {
std::cout << aRay.mStart(0) << ' ' << aRay.mStart(1) << ' ' << aRay.mStart(2) << ' '
          << aRay.mDirection(0) << ' ' << aRay.mDirection(1) << ' ' << aRay.mDirection(2) << '\n';
    throw 0;
  }
}

bool Medium::hits(Ray const& aRay) {
  try {
    auto hit = mSolver.solve4x(aRay.mStart, aRay.mDirection, mObject.getX());
    return hit.mValid && mObject.hasPixel(hit.mValue);
  }
  catch(...) {
std::cout << aRay.mStart(0) << ' ' << aRay.mStart(1) << ' ' << aRay.mStart(2) << ' '
          << aRay.mDirection(0) << ' ' << aRay.mDirection(1) << ' ' << aRay.mDirection(2) << '\n';
    return false;
  }
}


Image::Image(Parameters const& aPara, Medium &aMedium)
  : mRestrictCpu(aPara.mRestrictCpu)
  , mBuffer(aPara.mResolution * aPara.mResolution)
  , mImage(aPara.mResolution, aPara.mResolution)
  , mSubSample(aPara.mSubsample)
  , mSsFactor(1.0 / aPara.mSubsample)
  , mPixelSize(csFilmSize / aPara.mResolution)
  , mCenter(0.0, aPara.mCamCenter, 0.0)
  , mNormal(::cos(aPara.mTilt * cgPi / 180.0), ::sin(aPara.mTilt * cgPi / 180.0), 0.0)
  , mInPlaneZ(0.0, 0.0, 1.0)
  , mInPlaneY(mNormal.cross(mInPlaneZ))
  , mPinhole(mCenter + aPara.mPinholeDist * mNormal)
  , mBiasZ((aPara.mResolution - 1.0) / 2.0)
  , mBiasY((aPara.mResolution - 1.0) / 2.0)
  , mBiasSub((aPara.mSubsample - 1.0) / 2.0)
  , mColorGrid(csColorGrid * aPara.mSubsample * aPara.mSubsample)
  , mColorMirror(csColorMirror * aPara.mSubsample * aPara.mSubsample)
  , mGridIndent(std::max(0.0, std::min(1.0, aPara.mGridIndent)))
  , mGridSpacing(aPara.mGridSpacing)
  , mMirrorAcross(aPara.mMirrorAcross)
  , mMedium(aMedium) {}

void Image::process(char const * const aNameSurf, char const * const aNameOut) {
  calculateLimits();
  int mirrorHeight = calculateMirrorHeight();
  calculateMirage(mirrorHeight);
  if(*aNameSurf != 0) {
    renderSurface(aNameSurf);
  }
  else {} // nothing to do
  mImage.write(aNameOut);
}

void Image::calculateLimits() {
  Ray ray;
  ray.mStart = mPinhole;
  ray.mDirection = getDirectionInXy(csLimitLow - csLimitDelta);
  auto lastHit = mMedium.hits(ray);
  bool was = false;
  double limitAnglePrev;
  for(auto angle = csLimitLow; angle <= csLimitHigh; angle += csLimitDelta) {
    ray.mDirection = getDirectionInXy(angle);
    auto thisHit = mMedium.hits(ray);
    if(lastHit != thisHit) {
      auto critical = binarySearch(angle - csLimitDelta, angle, csLimitEpsilon, [this, &ray](auto const search){
        ray.mDirection = getDirectionInXy(search);
        return mMedium.hits(ray);
      });
      critical += (thisHit ? csLimitEpsilon : 0.0);
      std::cout << (thisHit ? "enter: " : "leave: ") << std::setprecision(10) << (critical * 180.0 / cgPi) << '\n';
      limitAnglePrev = mLimitAngleTop;
      mLimitAngleTop = critical * csLimitAngleBoost;
      if(!was) {
        mLimitAngleBottom = critical * csLimitAngleBoost;
        was = true;
      }
      else {} // nothing to do
    }
    else {} // nothing to do
    lastHit = thisHit;
  }
  auto angleY = (limitAnglePrev + mLimitAngleTop) / 2.0;
  ray.mDirection = getDirectionYz(angleY, csLimitLow - csLimitDelta);
  mLimitAngleDeep = binarySearch(csLimitLow, 0.0, csLimitEpsilon, [this, &ray, angleY](auto const search){
    ray.mDirection = getDirectionYz(angleY, search);
    return mMedium.hits(ray);
  });
  mLimitAngleDeep *= csLimitAngleBoost;
  mLimitAngleShallow = -mLimitAngleDeep;

  int y = mImage.get_height() / 2;
  int z;
  for(z = 0; z < mImage.get_width() / 2; ++z) {
    Vertex subpixel = mCenter + mPixelSize * (
      (z - mBiasZ) * mInPlaneZ +
      (y - mBiasY) * mInPlaneY);
    ray.mDirection = (mPinhole - subpixel).normalized();
    auto angleZ = -std::atan(ray.mDirection(2) / ray.mDirection(0));
    if(angleZ > mLimitAngleDeep) {
      break;
    }
    else {} // nothing to do
  }
  mLimitPixelDeep = z;
  mLimitPixelShallow = mImage.get_width() - z - 1;

  z = mImage.get_width() / 2;
  for(y = 0; y < mImage.get_height(); ++y) {
    Vertex subpixel = mCenter + mPixelSize * (
      (z - mBiasZ) * mInPlaneZ +
      (y - mBiasY) * mInPlaneY);
    ray.mDirection = (mPinhole - subpixel).normalized();
    auto angleY = std::atan(ray.mDirection(1) / ray.mDirection(0));
    if(angleY > mLimitAngleBottom) {
      break;
    }
    else {} // nothing to do
  }
  mLimitPixelBottom = y;  // surface rendering goes from 0 to this
}

int Image::calculateMirrorHeight() {
  int result = -1;
  double minHit = std::numeric_limits<double>::max();
  int z = mImage.get_width() / 2;
  Ray ray;
  ray.mStart = mPinhole;
  for(int y = 0; y < mImage.get_height(); ++y) {
    Vertex subpixel = mCenter + mPixelSize * (
      (z - mBiasZ) * mInPlaneZ +
      (y - mBiasY) * mInPlaneY);
    ray.mDirection = (mPinhole - subpixel).normalized();
    auto hit = mMedium.getHit(ray);
    if(hit.mValid && hit.mValue(1) < minHit) {
      minHit = hit.mValue(1);
      result = y;
    }
    else {} // nothing to do
  }
  return result;
}

void Image::calculateMirage(int const aMirrorHeight) {
  uint32_t nCpus = std::thread::hardware_concurrency();
  nCpus -= (nCpus <= mRestrictCpu ? nCpus - 1u : mRestrictCpu);
  std::vector<std::thread> threads(nCpus);
  for (uint32_t i = 0u; i < nCpus; ++i) {
    threads[i] = std::thread([this, aMirrorHeight, nCpus, i] {
      Ray ray;
      ray.mStart = mPinhole;
      Medium localMedium(mMedium);
      auto yBegin = i * mImage.get_height() / nCpus;
      auto yEnd = (i + 1u) * mImage.get_height() / nCpus;
      for(int y = yBegin; y < yEnd; ++y) {
        for(int z = 0; z < mImage.get_width(); ++z) {
          double sum = 0.0;
          for(uint32_t i = 0; i < mSubSample; ++i) {
            for(uint32_t j = 0; j < mSubSample; ++j) {
              Vertex subpixel = mCenter + mPixelSize * (
                    (z - mBiasZ + mSsFactor * (i - mBiasSub)) * mInPlaneZ +
                    (y - mBiasY + mSsFactor * (j - mBiasSub)) * mInPlaneY);
              ray.mDirection = (mPinhole - subpixel).normalized();
              auto angleY = std::atan(ray.mDirection(1) / ray.mDirection(0));
              auto angleZ = std::atan(ray.mDirection(2) / ray.mDirection(0));
              if(angleY > mLimitAngleBottom && angleY < mLimitAngleTop && angleZ > mLimitAngleDeep && angleZ < mLimitAngleShallow) {
                sum += localMedium.trace(ray);
              }
              else {} // nothing to do
            }
          }
          if(y == aMirrorHeight && mMirrorAcross) {
             sum = mColorMirror;
          }
          else {} // nothing to do
          if(z < mLimitPixelDeep * mGridIndent || z > mImage.get_width() - mLimitPixelDeep * mGridIndent) {
            if(y == aMirrorHeight) {
               sum = mColorMirror;
            }
            else if(mGridSpacing > 1u && y % mGridSpacing == 0) {
              sum = mColorGrid;
            }
            else {} // nothing to do
          }
          else {} // nothing to do
          mBuffer[(mImage.get_width() - z - 1u) + mImage.get_width() * (mImage.get_height() - y - 1u)] = ::round(sum / static_cast<double>(mSubSample * mSubSample));
        }
      }
    });
  }
  for (auto& t : threads) {
    t.join();
  }
  for(int y = 0; y < mImage.get_height(); ++y) {
    for(int z = 0; z < mImage.get_width(); ++z) {
      mImage.set_pixel(z, y, mBuffer[y * mImage.get_width() + z]);
    }
  }
}

void Image::renderSurface(char const * const aNameSurf) {
  auto ssFactor = 1.0 / csSurfSubsample;
  Vector normal(1.0, 0.0, 0.0);
  Vector inPlaneY(mNormal.cross(mInPlaneZ));
  Vertex pinhole(csSurfPinholeDist * normal);
  double biasSub = ((csSurfSubsample - 1.0) / 2.0);

  Plane plane = Plane::createFrom2vectors1point(Vector(0.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0), Vertex(csSurfaceDistance, 0.0, 0.0));
  Ray ray;
  ray.mStart = pinhole;

  Vertex pixel = mPixelSize * ((mLimitPixelDeep - mBiasZ) * mInPlaneZ + (mLimitPixelBottom - mBiasY) * inPlaneY);
  ray.mDirection = (pinhole - pixel).normalized();
  auto intersection = plane.intersect(ray);
  auto upperCorner = intersection.mPoint;

  pixel = mPixelSize * ((mLimitPixelShallow - mBiasZ) * mInPlaneZ - mBiasY * inPlaneY);
  ray.mDirection = (pinhole - pixel).normalized();
  intersection = plane.intersect(ray);
  auto lowerCorner = intersection.mPoint;

  Object surface(aNameSurf, upperCorner, lowerCorner);

  for(int y = 0; y < mLimitPixelBottom; ++y) {
    for(int z = mLimitPixelDeep; z < mLimitPixelShallow; ++z) {
      double sum = 0.0;
      for(uint32_t i = 0; i < csSurfSubsample; ++i) {
        for(uint32_t j = 0; j < csSurfSubsample; ++j) {
          Vertex subpixel = mPixelSize * (
                (z - mBiasZ + ssFactor * (i - biasSub)) * mInPlaneZ +
                (y - mBiasY + ssFactor * (j - biasSub)) * inPlaneY);
          ray.mDirection = (pinhole - subpixel).normalized();
          auto intersection = plane.intersect(ray);
          sum += surface.getPixel(intersection.mPoint);
        }
      }
      mImage.set_pixel(mImage.get_width() - z - 1u, mImage.get_height() - y - 1u, ::round(sum / static_cast<double>(csSurfSubsample * csSurfSubsample)));
    }
  }
}
