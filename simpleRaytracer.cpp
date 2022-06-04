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
std::cout << mDy << ' ' << mDz << '\n';
  std::cout << "minY: " << mMinY << '\n';
  std::cout << "maxY: " << mMaxY << '\n';
  std::cout << "minZ: " << mMinZ << '\n';
  std::cout << "maxZ: " << mMaxZ << '\n';
}

Object::Object(char const * const aName, Vertex const& aUpperCorner, Vertex const& aLowerCorner)
  : mImage(aName)
  , mDy((aUpperCorner(2) - aLowerCorner(2)) / mImage.get_width())
  , mDz((aUpperCorner(2) - aLowerCorner(2)) / mImage.get_width())
  , mMinY(aUpperCorner(1))
  , mMaxY(aUpperCorner(1) + ((aUpperCorner(2) - aLowerCorner(2)) * mImage.get_height()) / mImage.get_width())
  , mMinZ(aLowerCorner(2))
  , mMaxZ(aUpperCorner(2))
  , mX(aLowerCorner(0)) {
std::cout << aUpperCorner(1) << '\n';
std::cout << aUpperCorner(2) << '\n';
std::cout << aLowerCorner(2) << '\n';
std::cout << mImage.get_height() << '\n';
std::cout << mImage.get_width() << '\n';
std::cout << mDy << ' ' << mDz << '\n';
  std::cout << "minY: " << mMinY << '\n';
  std::cout << "maxY: " << mMaxY << '\n';
  std::cout << "minZ: " << mMinZ << '\n';
  std::cout << "maxZ: " << mMaxZ << '\n';
}

bool Object::hasPixel(Vertex const &aHit) const {
  return aHit(1) > mMinY && aHit(1)  < mMaxY && aHit(2) > mMinZ && aHit(2) < mMaxZ;
}

uint8_t Object::getPixel(Vertex const &aHit) const {
  uint8_t result = 0u;
  int32_t x = static_cast<int32_t>(::round((aHit(2) - mMinZ) / mDz));
  int32_t y = static_cast<int32_t>(::round(mImage.get_height() - (aHit(1) - mMinY) / mDy - 1u));
  if(x >= 0 && y >= 0 && x < mImage.get_width() && y < mImage.get_height()) {
    result = mImage.get_pixel(x, y);
//std::cout << "h " << aHit(2) << ' ' << aHit(1) << "  /  " << x << ' ' << y << "  -  " << result << '\n';
  }
  else {
// std::cout << "M " << aHit(2) << ' ' << aHit(1) <<  "  /  " << x << ' ' << y <<'\n';
} // Nothing to do
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
  , mPalette(256)
  , mResolutionX(aPara.mResolutionX)
  , mBorderFactor(aPara.mBorderFactor)
  , mSubSample(aPara.mSubsample)
  , mSsFactor(1.0 / aPara.mSubsample)
  , mCenter(0.0, aPara.mCamCenter, 0.0)
  , mNormal(::cos(aPara.mTilt * cgPi / 180.0), ::sin(aPara.mTilt * cgPi / 180.0), 0.0)
  , mInPlaneZ(0.0, 0.0, 1.0)
  , mInPlaneY(mNormal.cross(mInPlaneZ))
  , mPinhole(mCenter + csSurfPinholeDist * mNormal)
  , mBiasSub((aPara.mSubsample - 1.0) / 2.0)
  , mMarkIndent(std::max(0.0, std::min(1.0, aPara.mMarkIndent)))
  , mMarkAcross(aPara.mMarkAcross)
  , mMedium(aMedium) {
  mPalette[csColorMirror] = png::color(255u, 0u, 0u);
  mPalette[csColorBase] = png::color(0u, 255u, 0u);
  for(uint32_t i = csColorBlack; i < mPalette.size(); ++i) {
    mPalette[i] = png::color(i, i, i);
  }
  mImage.set_palette(mPalette);
}

void Image::process(char const * const aNameSurf, char const * const aNameOut) {
  calculateAngleLimits(Eikonal::Temperature::cAmbient);
  calculateAngleLimits(Eikonal::Temperature::cBase);
  calculateAngleLimits(Eikonal::Temperature::cMinimum);
  calculateAngleLimits(Eikonal::Temperature::cMaximum);
  calculateBiases(*aNameSurf != 0);
std::cout << "l angleTop:     " <<  *mLimitAngleTop << '\n';
std::cout << "l angleBottom:  " <<  *mLimitAngleBottom << '\n';
std::cout << "l angleDeep:     " <<  *mLimitAngleDeep << '\n';
std::cout << "l angleShallow:  " << *mLimitAngleShallow << '\n';
  mLimitAngleTop.reset();
  mLimitAngleBottom.reset();
  calculateAngleLimits(Eikonal::Temperature::cAmbient);
  mLimitPixelBottom     = calculatePixelLimitY(*mLimitAngleBottom);
std::cout << "a angleBottom:  " <<  *mLimitAngleBottom << '\n';
std::cout << "bottom:     " << mLimitPixelBottom << '\n';
  mLimitAngleTop.reset();
  mLimitAngleBottom.reset();
  calculateAngleLimits(Eikonal::Temperature::cBase);
std::cout << "b angleTop:  " <<  *mLimitAngleTop << '\n';
std::cout << "b angleBottom:  " <<  *mLimitAngleBottom << '\n';
std::cout << "b angleBottomSurf:  " <<  mLimitAngleBottomSurf << '\n';
  mLimitPixelBaseTop        = calculatePixelLimitY(*mLimitAngleTop);
  mLimitPixelBaseBottom     = calculatePixelLimitY(*mLimitAngleBottom);
  mLimitPixelBaseBottomSurf = calculatePixelLimitY(mLimitAngleBottomSurf);
  calculateAngleLimits(Eikonal::Temperature::cMinimum);
  calculateAngleLimits(Eikonal::Temperature::cMaximum);
  calculateAngleLimits(Eikonal::Temperature::cAmbient);
std::cout << "l angleTop:     " <<  *mLimitAngleTop << '\n';
std::cout << "l angleBottom:  " <<  *mLimitAngleBottom << '\n';
std::cout << "l angleDeep:     " <<  *mLimitAngleDeep << '\n';
std::cout << "l angleShallow:  " << *mLimitAngleShallow << '\n';
  mLimitPixelDeep       = calculatePixelLimitZ(*mLimitAngleDeep);
  mLimitPixelShallow    = calculatePixelLimitZ(*mLimitAngleShallow);
std::cout << "baseTop:        " << mLimitPixelBaseTop << '\n';
std::cout << "baseBottom:     " << mLimitPixelBaseBottom << '\n';
std::cout << "baseBottomSurf: " << mLimitPixelBaseBottomSurf << '\n';
std::cout << "deep:           " << mLimitPixelDeep << '\n';
std::cout << "shallow:        " << mLimitPixelShallow << '\n';
std::cout << "resX:           " << mImage.get_width() << '\n';
std::cout << "resY:           " << mImage.get_height() << '\n';
std::cout << "biasZ:          " << mBiasZ << '\n';
std::cout << "biasY:          " << mBiasY << '\n';
  if(*aNameSurf != 0) {
    renderSurface(aNameSurf);
  }
  else {} // nothing to do
std::cout << "errfdddddddddddddddddddddddddd\n";
  int mirrorHeight = calculateMirrorHeight();
  calculateMirage();
  drawMarks(mirrorHeight);
  mImage.write(aNameOut);
}

void Image::calculateAngleLimits(Eikonal::Temperature const aWhich) {
  mMedium.setWaterTempAmb(aWhich);
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
      limitAnglePrev = mLimitAngleTop.value_or(0.0);
      auto tmp = critical * csLimitAngleBoost;
      mLimitAngleTop = (mLimitAngleTop ? std::max(*mLimitAngleTop, tmp) : tmp);
      if(!was) {
        mLimitAngleBottom = (mLimitAngleBottom ? std::min(*mLimitAngleBottom, tmp) : tmp);
        was = true;
      }
      else {} // nothing to do
    }
    else {} // nothing to do
    lastHit = thisHit;
  }
  auto angleY = (limitAnglePrev + *mLimitAngleTop) / 2.0;
  ray.mDirection = getDirectionYz(angleY, csLimitLow - csLimitDelta);
  auto tmp = binarySearch(csLimitLow, 0.0, csLimitEpsilon, [this, &ray, angleY](auto const search){
    ray.mDirection = getDirectionYz(angleY, search);
    return mMedium.hits(ray);
  });
  tmp *= csLimitAngleBoost;
  mLimitAngleDeep = (mLimitAngleDeep ? std::min(*mLimitAngleDeep, tmp) : tmp);
  mLimitAngleShallow = -*mLimitAngleDeep;
}

void Image::calculateBiases(bool const aRenderSurface) {
  auto film = Plane::createFrom2vectors1point(mInPlaneY, mInPlaneZ, mCenter);
  auto angleDiff    = *mLimitAngleTop + *mLimitAngleShallow - *mLimitAngleBottom - *mLimitAngleDeep;

  auto angleTop               = *mLimitAngleTop     + angleDiff * mBorderFactor;
  auto angleBottom            = *mLimitAngleBottom  - angleDiff * mBorderFactor * (aRenderSurface ? csRenderSurfaceFactor : 1.0);
       mLimitAngleBottomSurf  = *mLimitAngleBottom  - angleDiff * mBorderFactor * (csRenderSurfaceFactor - 1.0);
  auto angleDeep              = *mLimitAngleDeep    - angleDiff * mBorderFactor;
  auto angleShallow           = *mLimitAngleShallow + angleDiff * mBorderFactor;

  auto limitTop     = film.intersect(mPinhole, -getDirectionInXy(angleTop)).mPoint;
  auto limitBottom  = film.intersect(mPinhole, -getDirectionInXy(angleBottom)).mPoint;
  auto limitDeep    = film.intersect(mPinhole, -getDirectionInXz(angleDeep)).mPoint;
  auto limitShallow = film.intersect(mPinhole, -getDirectionInXz(angleShallow)).mPoint;

  auto height = (limitTop - limitBottom).norm();
  auto width  = (limitDeep - limitShallow).norm();
  auto resolutionY = static_cast<uint32_t>(std::round(height / width * mResolutionX));
  resolutionY += resolutionY % 2u; // To allow an MPEG converter run on the output.

  mBiasZ = (mResolutionX - 1.0) * (mCenter - limitDeep).norm() / width;
  mBiasY = (resolutionY - 1.0) * (mCenter - limitBottom).norm() / height;
  mPixelSize = (width / mResolutionX + height / resolutionY) / 2.0;

  mBuffer.reserve(mResolutionX * resolutionY);
  mBuffer.insert(mBuffer.begin(), mResolutionX * resolutionY, csColorVoid);
  mImage.resize(mResolutionX, resolutionY);
}

int Image::calculatePixelLimitZ(double const aAngle) {
  Ray ray;
  ray.mStart = mPinhole;
  int y = mImage.get_height() / 2;
  int result;
  double minDist = std::numeric_limits<double>::max();
  for(int z = 0; z < mImage.get_width(); ++z) {
    Vertex subpixel = mCenter + mPixelSize * (
      (z - mBiasZ) * mInPlaneZ +
      (y - mBiasY) * mInPlaneY);
    ray.mDirection = (mPinhole - subpixel).normalized();
    auto angleZ = -std::atan(ray.mDirection(2) / ray.mDirection(0));
    auto now = std::abs(angleZ - aAngle);
    if(now < minDist) {
      minDist = now;
      result = z;
    }
    else {} // nothing to do
  }
  return result;
}

int Image::calculatePixelLimitY(double const aAngle) {
  Ray ray;
  ray.mStart = mPinhole;
  int z = mImage.get_width() / 2;
  int result;
  double minDist = std::numeric_limits<double>::max();
  for(int y = 0; y < mImage.get_height(); ++y) {
    Vertex subpixel = mCenter + mPixelSize * (
      (z - mBiasZ) * mInPlaneZ +
      (y - mBiasY) * mInPlaneY);
    ray.mDirection = (mPinhole - subpixel).normalized();
    auto angleY = std::atan(ray.mDirection(1) / ray.mDirection(0));
    auto now = std::abs(angleY - aAngle);
    if(now < minDist) {
      minDist = now;
      result = y;
    }
    else {} // nothing to do
  }
  return result;
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

void Image::renderSurface(char const * const aNameSurf) {
  auto ssFactor = 1.0 / csSurfSubsample;
  double biasSub = ((csSurfSubsample - 1.0) / 2.0);

  Plane plane = Plane::createFrom2vectors1point(Vector(0.0, 1.0, 0.0), Vector(0.0, 0.0, 1.0), Vertex(csSurfaceDistance, 0.0, 0.0));
  Ray ray;
  ray.mStart = mPinhole;

std::cout << "bottom:     " << mLimitPixelBottom << '\n';
std::cout << "baseBottomSurf: " << mLimitPixelBaseBottomSurf << '\n';
  Vertex pixel = mCenter + mPixelSize * ((mLimitPixelDeep - mBiasZ) * mInPlaneZ + (mLimitPixelBaseBottomSurf - mBiasY) * mInPlaneY);
  ray.mDirection = (mPinhole - pixel).normalized();
  auto intersection = plane.intersect(ray);
  auto upperCorner = intersection.mPoint;

  pixel = mCenter + mPixelSize * ((mLimitPixelShallow - mBiasZ) * mInPlaneZ - mBiasY * mInPlaneY);
  ray.mDirection = (mPinhole - pixel).normalized();
  intersection = plane.intersect(ray);
  auto lowerCorner = intersection.mPoint;

  Object surface(aNameSurf, upperCorner, lowerCorner);

  for(int y = mLimitPixelBaseBottomSurf; y <= mLimitPixelBottom; ++y) {
    for(int z = mLimitPixelDeep; z < mLimitPixelShallow; ++z) {
      double sum = 0.0;
      for(uint32_t i = 0; i < csSurfSubsample; ++i) {
        for(uint32_t j = 0; j < csSurfSubsample; ++j) {
          Vertex subpixel = mCenter + mPixelSize * (
                (z - mBiasZ + ssFactor * (i - biasSub)) * mInPlaneZ +
                (y - mBiasY + ssFactor * (j - biasSub)) * mInPlaneY);
          ray.mDirection = (mPinhole - subpixel).normalized();
          auto intersection = plane.intersect(ray);
          sum += surface.getPixel(intersection.mPoint);
        }
      }
      uint8_t color = std::max(csColorBlack, static_cast<uint8_t>(::round(sum / static_cast<double>(csSurfSubsample * csSurfSubsample))));
      mImage.set_pixel(mImage.get_width() - z - 1u, mImage.get_height() - y - 1u, color);
    }
  }
}

void Image::calculateMirage() {
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
          bool was = false;
          for(uint32_t i = 0; i < mSubSample; ++i) {
            for(uint32_t j = 0; j < mSubSample; ++j) {
              Vertex subpixel = mCenter + mPixelSize * (
                    (z - mBiasZ + mSsFactor * (i - mBiasSub)) * mInPlaneZ +
                    (y - mBiasY + mSsFactor * (j - mBiasSub)) * mInPlaneY);
              ray.mDirection = (mPinhole - subpixel).normalized();
              auto angleY = std::atan(ray.mDirection(1) / ray.mDirection(0));
              auto angleZ = std::atan(ray.mDirection(2) / ray.mDirection(0));
              if(angleY > mLimitAngleBottom && angleY < mLimitAngleTop && angleZ > mLimitAngleDeep && angleZ < mLimitAngleShallow) {
                was = true;
                sum += localMedium.trace(ray);
              }
              else {} // nothing to do
            }
          }
          uint8_t color;
          if(was) {
            color = std::max(csColorBlack, static_cast<uint8_t>(::round(sum / static_cast<double>(mSubSample * mSubSample))));
          }
          else {
            color = csColorVoid;
          }
          mBuffer[(mImage.get_width() - z - 1u) + mImage.get_width() * (mImage.get_height() - y - 1u)] = color;
        }
      }
    });
  }
  for (auto& t : threads) {
    t.join();
  }
  for(int y = 0; y < mImage.get_height(); ++y) {
    for(int z = 0; z < mImage.get_width(); ++z) {
      auto color = mBuffer[y * mImage.get_width() + z];
      if(color != csColorVoid) {
        mImage.set_pixel(z, y, color);
      }
      else {} // nothing to do
    }
  }
}

void Image::drawMarks(int const aMirrorHeight) {
  auto dashLength = std::max(static_cast<int>(mImage.get_width() / csDashCount), 2);
  auto dashLimit  = dashLength / 2;
  for(int z = 0; z < mImage.get_width(); ++z) {
    if(mMarkAcross || z < mLimitPixelDeep * mMarkIndent || z > mImage.get_width() - mLimitPixelDeep * mMarkIndent) {
      if(aMirrorHeight >= 0 && aMirrorHeight < mImage.get_width() && (z % dashLength < dashLimit)) {
        mImage.set_pixel(z, mImage.get_height() - 1 - aMirrorHeight, csColorMirror);
      }
      else {} // nothing to do
      if(mLimitPixelBaseTop >= 0 && mLimitPixelBaseTop < mImage.get_width() && (z % dashLength >= dashLimit)) {
        mImage.set_pixel(z, mImage.get_height() - 1 - mLimitPixelBaseTop, csColorBase);
      }
      else {} // nothing to do
      if(mLimitPixelBaseBottom >= 0 && mLimitPixelBaseBottom < mImage.get_width() && (z % dashLength >= dashLimit)) {
        mImage.set_pixel(z, mImage.get_height() - 1 - mLimitPixelBaseBottom, csColorBase);
      }
      else {} // nothing to do
    }
    else {} // nothing to do
  }
}
