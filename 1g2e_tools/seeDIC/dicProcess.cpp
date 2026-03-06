#include "dicProcess.hpp"

void DICResult::read(const char* path, int num) {
  char fname[256];
  sprintf(fname, "%s/dic_out_%d.txt", path, num);
  ifstream file(fname);
  if (!file) {
    cerr << "A problem occurred when opening the file " << fname << endl;
  }
  grain G;
  int ng = 0;
  grains.clear();

  file >> ng;
  for (int i = 0; i < ng; i++) {
    file >> G.refcoord_xpix >> G.refcoord_ypix >> G.refrot >> G.radius_pix >> G.dx >> G.dy >> G.drot >> G.upix >>
        G.vpix >> G.rot_inc >> G.NCC >> G.NCC_rescue >> G.NCC_subpix;
    grains.push_back(G);
  }
}
