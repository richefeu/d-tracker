#ifndef DICPROCESS_HPP_6BB571C4
#define DICPROCESS_HPP_6BB571C4

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

struct AABB {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  // This corresponds to the size we get with PhaseOne IQ80.
  // But with the fit these values will be re-computed
  AABB() : xmin(0.0), xmax(10320.0), ymin(0.0), ymax(7752.0) { }
};


struct grain {
  double refcoord_xpix;
  double refcoord_ypix;
  double refrot;
  double radius_pix;          
  double dx;  
  double dy; 
  double drot;
  double upix;
  double vpix;
  double rot_inc;
  double NCC;
  double NCC_rescue;
  double NCC_subpix;
};

/// A configuration is the state of a packing at a given time.
/// It corresponds to a photographs for experimental data or
/// a configuration at a saved time for data from simulation.
struct DICResult {
  vector<grain> grains;

  DICResult() {}
  void read(const char* path, int num);
};

struct Commands {
  
  Commands() {}
  void read(const char* filename);
};

#endif /* end of include guard: DICPROCESS_HPP_6BB571C4 */