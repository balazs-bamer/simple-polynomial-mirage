#ifndef SHEPARDINTERPOLATION
#define SHEPARDINTERPOLATION

#include <array>
#include <limits>
#include <memory>
#include <variant>
#include <Eigen/Dense>


template<typename tCoordinate, typename tPayload, size_t tInPlace>
class Shepard3d final {
  using Location = Eigen::Array<tCoordinate, 3, 1>;

private:
  class Node final {
  private:
    using Children = std::array<std::unique_ptr<Node>, 8>;
    using Payload  = std::array<tPayload, tInPlace>;

    uint32_t                        mCountTotal;
    uint32_t                        mCountHere;
    std::variant<Children, Payload> mContents;
    Location                        mSize;
  };

  Location                             mBoundsMin;
  Location                             mBoundsMax;
  std::array<std::unique_ptr<Node>, 8> mRoots;

public:
  Shepard3d(std::vector<Location> const &aLocation, std::vector<tPayload> const &aData);

  tPayload interpolate(tCoordinate const aX, tCoordinate const aY, tCoordinate const aZ, uint32_t const aSamplesToConsider) const;
  tPayload interpolate(Location const& aLocation, uint32_t const aSamplesToConsider) const { return interpolate(aLocation[0], aLocation[1], aLocation[2], aSamplesToConsider); }
};

template<typename tCoordinate, typename tPayload, size_t tInPlace>
Shepard3d<tCoordinate, tPayload, tInPlace>::Shepard3d(std::vector<Location> const &aLocation, std::vector<tPayload> const &aData)
 : mBoundsMin {  std::numeric_limits<tCoordinate>::max(),  std::numeric_limits<tCoordinate>::max(),  std::numeric_limits<tCoordinate>::max()}
 , mBoundsMax { -std::numeric_limits<tCoordinate>::max(), -std::numeric_limits<tCoordinate>::max(), -std::numeric_limits<tCoordinate>::max()} {
  if(aLocation.size() == aData.size() && aData.size() > 0u) {
    for(size_t i = 0u; i < aLocation.size(); ++i) {
      for(size_t j = 0u; j < 3u; ++j) {
        mBoundsMin[j] = std::min(mBoundsMin[j], aLocation[i][j]);
        mBoundsMax[j] = std::max(mBoundsMax[j], aLocation[i][j]);
      }
    }
  }
  else {} // Nothing to do
}

#endif // SHEPARD_H
