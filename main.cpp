#include "Angle2apparentMirrorDepth.h"
#include <iostream>

int main(int argc, char **argv) {
/*  Angle2apparentMirrorDepth t32(20.0f, 12.0f);
  Angle2apparentMirrorDepth t40(20.0f, 20.0f);
  Angle2apparentMirrorDepth t48(20.0f, 28.0f);
  Angle2apparentMirrorDepth t56(20.0f, 36.0f);*/
  Angle2apparentMirrorDepth t64(20.0f, 44.0f);
  std::array<float, 12> h = {{0.0005f, 0.01f, 0.03, 0.07f, 0.13f, 0.25f, 0.38f, 0.55f, 0.65f, 0.77f, 0.89f, 1.0f}};
  for(auto f : h) {
    std::cout << "h " << f;
/*    std::cout << " t32 " << t32.getTempAtHeight(f);
    std::cout << " t40 " << t40.getTempAtHeight(f);
    std::cout << " t48 " << t48.getTempAtHeight(f);
    std::cout << " t56 " << t56.getTempAtHeight(f);*/
    std::cout << " t64 " << t64.getTempAtHeight(f);
    std::cout << '\n';
  }
  return 0;
}
