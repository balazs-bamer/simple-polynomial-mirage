#ifndef SHEPARDINTERPOLATION
#define SHEPARDINTERPOLATION

#include <array>
#include <limits>
#include <memory>
#include <variant>
#include <Eigen/Dense>


template<typename tCoordinate, typename tPayload, size_t tInPlace>
class Shepard3d final {
public:
  using Location = Eigen::Array<tCoordinate, 3, 1>;
  struct Data {
    Location mLocation;
    tPayload mPayload;
  };

private:
  struct Node final {
  public:
    using Children = std::array<std::unique_ptr<Node>, 8>;
    using Payload  = std::array<Data, tInPlace>;

    uint32_t                        mCountTotal;
    uint32_t                        mCountHere;
    std::variant<Children, Payload> mContents;
    Location                        mCenter;
    Location                        mSize;

    Node(Location const &aCenter, Location const& aSize) : mCountTotal(0u), mCountHere(0u), mCenter(aCenter), mSize(aSize) {}
  };

  std::array<std::unique_ptr<Node>, 8> mRoots;
  Location                             mBoundsMin;
  Location                             mBoundsMax;

public:
  Shepard3d(std::vector<Data> const &aData, uint32_t const aSamplesToConsider);

  tPayload interpolate(tCoordinate const aX, tCoordinate const aY, tCoordinate const aZ) const;
  tPayload interpolate(Location const& aLocation) const { return interpolate(aLocation[0], aLocation[1], aLocation[2]); }

private:
  void buildTree(size_t const aWhich, Location const aCenter, Location const aBoundsMax, std::vector<Data> const& aData);
};

template<typename tCoordinate, typename tPayload, size_t tInPlace>
Shepard3d<tCoordinate, tPayload, tInPlace>::Shepard3d(std::vector<Data> const &aData, uint32_t const aSamplesToConsider)
  : mBoundsMin {  std::numeric_limits<tCoordinate>::max(),  std::numeric_limits<tCoordinate>::max(),  std::numeric_limits<tCoordinate>::max()}
  , mBoundsMax { -std::numeric_limits<tCoordinate>::max(), -std::numeric_limits<tCoordinate>::max(), -std::numeric_limits<tCoordinate>::max()} {
  for(auto const &item : aData) {
    for(size_t j = 0u; j < 3u; ++j) {
      mBoundsMin[j] = std::min(mBoundsMin[j], item.mLocation[j]);
      mBoundsMax[j] = std::max(mBoundsMax[j], item.mLocation[j]);
    }
  }

  buildTree(0u, (mBoundsMin + mBoundsMax) / 2u, mBoundsMax - mBoundsMin, &aData);
}

template<typename tCoordinate, typename tPayload, size_t tInPlace>
void Shepard3d<tCoordinate, tPayload, tInPlace>::buildTree(size_t const aWhich, Location const aCenter, Location const aSize, std::vector<Data> const& aData) {
  mRoots[aWhich] = std::move(std::make_unique<Node>(aCenter, aSize));
  auto root = mRoots[aWhich].get();
  for(auto const &item : aData) {
    auto branch = root;
    while(true) {
      ++branch->mCountTotal;
      if(branch->mCountHere < tInPlace) {
        if(branch->mCountHere == 0u) {
          branch->mContents = Node::Payload();
        }
        else {} // nothing to do
        std::get<Node::Payload>(branch->mContents)[branch->mCountHere] = item;
        ++branch->mCountHere;
        break;
      }
      else {

      }
    }
  }
}

#endif // SHEPARD_H
