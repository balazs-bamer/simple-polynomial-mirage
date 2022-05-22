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
  std::cout << mMinZ << ' ' << mMaxZ << ' ' << mMinY << ' ' << mMaxY << '\n';
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


bool eq(double a, double b) { return std::abs(a-b)<1e-7; }

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
  , mGridColor(static_cast<double>(std::max(0u, std::min(255u, aPara.mGridColor))) * aPara.mSubsample * aPara.mSubsample)
  , mGridIndent(std::max(0.0, std::min(1.0, aPara.mGridIndent)))
  , mGridSpacing(aPara.mGridSpacing)
  , mMedium(aMedium) {}

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
  mLimitPixelDeep = static_cast<int>(z * mGridIndent);
  mLimitPixelShallow = mImage.get_width() - mLimitPixelDeep;
}

void Image::process(char const * const aName) {
  uint32_t nCpus = std::thread::hardware_concurrency();
  nCpus -= (nCpus <= mRestrictCpu ? nCpus - 1u : mRestrictCpu);
  std::vector<std::thread> threads(nCpus);
  for (uint32_t i = 0u; i < nCpus; ++i) {
    threads[i] = std::thread([this, nCpus, i] {
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
          if(mGridSpacing > 1u && y % mGridSpacing == 0 && (z < mLimitPixelDeep || z > mLimitPixelShallow)) {
            sum = mGridColor;
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
  mImage.write(aName);
}
