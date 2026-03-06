#ifndef CONFPROCESS_HPP_6BB571C4
#define CONFPROCESS_HPP_6BB571C4

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

using namespace std;

const int format_x_y_R_rot = 0;
const int format_x_y_R_rot_NCC_NCC = 1;
const int format_simuDEM = 2;

struct grain {
  // instantaneous data
  double x, y, R;       // position and radius
  double rot;          // orientation (angle in radian, anticlockwise, 0 corresponds to horizontal)
  double zncc1, zncc2;  // ZNCCs before and after subpixel tracking

  // data computed from reference to current configuration
  double dx, dy, drot;     // displacement
  double ddx, ddy, ddrot;  // displacement fluctuation
  double exx, eyy, exy;    // mean local strain
};

struct AABB {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
};

//   0 +-------+ 1
//    /       /
// 3 +-------+ 2
struct box {
  double x[4], y[4];

  /// Returns the Axis Aligned Bounding Box (AABB) that surrounds the current box
  AABB getAABB() {
    AABB aabb;
    aabb.xmin = x[0];
    aabb.xmax = x[0];
    aabb.ymin = y[0];
    aabb.ymax = y[0];
    for (size_t i = 1; i < 4; i++) {
      if (aabb.xmin > x[i]) aabb.xmin = x[i];
      if (aabb.xmax < x[i]) aabb.xmax = x[i];
      if (aabb.ymin > y[i]) aabb.ymin = y[i];
      if (aabb.ymax < y[i]) aabb.ymax = y[i];
    }
    return aabb;
  }
};

struct strain {
  double Gxx, Gyy, Gxy, Gyx, Exx, Eyy, Exy, Erot, Evol;
};

/// The neighbor struct holds the id-numbers of the bodies that are close.
/// The two bodies are not necessary in contact (but they can). The contact list
/// is thus a subset of the neighbor list.
/// Each neighbor holds also the branch length, and if there is a contact
/// it holds the position vector (x,y), the normal vector (nx,ny), fn and ft.
struct neighbor {
  size_t i, j;
  double x, y;
  double nx, ny;
  double fn, ft;
  double branch;

  neighbor() : fn(0.0), ft(0.0) {}  // Ctor
  neighbor(const neighbor& N) {     // Copy Ctor
    i = N.i;
    j = N.j;
    nx = N.nx;
    ny = N.ny;
    fn = N.fn;
    ft = N.ft;
    branch = N.branch;
  }
};

/// A configuration is the state of a packing at a given time.
/// It corresponds to a photographs for experimental data or
/// a configuration at a saved time for data from simulation.
struct Configuration {
  vector<grain> grains;
  vector<neighbor> neighbors;
  box Box;
  strain Strain;
  double scaleFactor;  // it is not clear if this is pixels/meter or meters/pixel

  // Ctor
  Configuration()
      : scaleFactor(1.0) {}

  void copy(const Configuration& C) {
    grains.clear();
    for (size_t i = 0; i < C.grains.size(); i++) grains.push_back(C.grains[i]);
    neighbors.clear();
    for (size_t i = 0; i < C.neighbors.size(); i++) neighbors.push_back(C.neighbors[i]);
    Box = C.Box;
    Strain = C.Strain;
    scaleFactor = C.scaleFactor;
  }

  void read(const char* path, int num, int format);
  void computeDisplacements(Configuration& RefConf);
  void computeFluctuations(Configuration& RefConf);
  void computeLocalMeanStrains(Configuration& RefConf);
  void updateNeighbors(double tol);
  double getRmin();
  double getRmax();
  double getRawOverlay(size_t ineighbor);
};

std::vector<std::pair<size_t, size_t> > buildNetworkCorrespondance(std::vector<neighbor>& A, std::vector<neighbor>& B);

#endif /* end of include guard: CONFPROCESS_HPP_6BB571C4 */