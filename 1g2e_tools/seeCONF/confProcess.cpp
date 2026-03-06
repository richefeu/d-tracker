#include "confProcess.hpp"

// ===========================================

int equals(double a, double b, double tolerance) {
  return (a == b) || ((a <= (b + tolerance)) && (a >= (b - tolerance)));
}

double cross2(double x0, double y0, double x1, double y1) { return x0 * y1 - y0 * x1; }

int in_range(double val, double range_min, double range_max, double tol) {
  return ((val + tol) >= range_min) && ((val - tol) <= range_max);
}

// Returns number of solutions found.  If there is one valid solution, it will be put in s and t
//   p2 --- p3
//   |      |
// t |  p   |
//   |      |
//   p0 --- p1
//      s
int inverseBilinear(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double x,
                    double y, double* sout, double* tout, double* s2out, double* t2out) {
  int t_valid, t2_valid;

  double a = cross2(x0 - x, y0 - y, x0 - x2, y0 - y2);
  double b1 = cross2(x0 - x, y0 - y, x1 - x3, y1 - y3);
  double b2 = cross2(x1 - x, y1 - y, x0 - x2, y0 - y2);
  double c = cross2(x1 - x, y1 - y, x1 - x3, y1 - y3);
  double b = 0.5 * (b1 + b2);

  double s, s2, t, t2;

  double am2bpc = a - 2 * b + c;
  // this is how many valid s values we have
  int num_valid_s = 0;

  if (equals(am2bpc, 0, 1e-10)) {
    if (equals(a - c, 0, 1e-10)) {
      // Looks like the input is a line
      // You could set s=0.5 and solve for t if you wanted to
      return 0;
    }
    s = a / (a - c);
    if (in_range(s, 0, 1, 1e-10)) num_valid_s = 1;
  } else {
    double sqrtbsqmac = sqrt(b * b - a * c);
    s = ((a - b) - sqrtbsqmac) / am2bpc;
    s2 = ((a - b) + sqrtbsqmac) / am2bpc;
    num_valid_s = 0;
    if (in_range(s, 0, 1, 1e-10)) {
      num_valid_s++;
      if (in_range(s2, 0, 1, 1e-10)) num_valid_s++;
    } else {
      if (in_range(s2, 0, 1, 1e-10)) {
        num_valid_s++;
        s = s2;
      }
    }
  }

  if (num_valid_s == 0) return 0;

  t_valid = 0;
  if (num_valid_s >= 1) {
    double tdenom_x = (1 - s) * (x0 - x2) + s * (x1 - x3);
    double tdenom_y = (1 - s) * (y0 - y2) + s * (y1 - y3);
    t_valid = 1;
    if (equals(tdenom_x, 0, 1e-10) && equals(tdenom_y, 0, 1e-10)) {
      t_valid = 0;
    } else {
      // Choose the more robust denominator
      if (fabs(tdenom_x) > fabs(tdenom_y)) {
        t = ((1 - s) * (x0 - x) + s * (x1 - x)) / (tdenom_x);
      } else {
        t = ((1 - s) * (y0 - y) + s * (y1 - y)) / (tdenom_y);
      }
      if (!in_range(t, 0, 1, 1e-10)) t_valid = 0;
    }
  }

  // Same thing for s2 and t2
  t2_valid = 0;
  if (num_valid_s == 2) {
    double tdenom_x = (1 - s2) * (x0 - x2) + s2 * (x1 - x3);
    double tdenom_y = (1 - s2) * (y0 - y2) + s2 * (y1 - y3);
    t2_valid = 1;
    if (equals(tdenom_x, 0, 1e-10) && equals(tdenom_y, 0, 1e-10)) {
      t2_valid = 0;
    } else {
      // Choose the more robust denominator
      if (fabs(tdenom_x) > fabs(tdenom_y)) {
        t2 = ((1 - s2) * (x0 - x) + s2 * (x1 - x)) / (tdenom_x);
      } else {
        t2 = ((1 - s2) * (y0 - y) + s2 * (y1 - y)) / (tdenom_y);
      }
      if (!in_range(t2, 0, 1, 1e-10)) t2_valid = 0;
    }
  }

  // Final cleanup
  if (t2_valid && !t_valid) {
    s = s2;
    t = t2;
    t_valid = t2_valid;
    t2_valid = 0;
  }

  // Output
  if (t_valid) {
    *sout = s;
    *tout = t;
  }

  if (t2_valid) {
    *s2out = s2;
    *t2out = t2;
  }

  return t_valid + t2_valid;
}

void bilinear(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3, double s,
              double t, double* x, double* y) {
  *x = t * (s * x3 + (1 - s) * x2) + (1 - t) * (s * x1 + (1 - s) * x0);
  *y = t * (s * y3 + (1 - s) * y2) + (1 - t) * (s * y1 + (1 - s) * y0);
}

// ===========================================

// Update the neighbor list.
void Configuration::updateNeighbors(double tolNeighbor) {
  // pack the current network
  vector<neighbor> packed;
  packed = neighbors;

  // update the new list
  neighbors.clear();
  neighbor C;
  for (size_t i = 0; i < grains.size(); i++) {
    for (size_t j = i + 1; j < grains.size(); j++) {
      double dx = grains[j].x - grains[i].x;
      double dy = grains[j].y - grains[i].y;
      double Ri = grains[i].R;
      double Rj = grains[j].R;
      double dd = dx * dx + dy * dy;
      double RR = Ri + Rj + tolNeighbor;
      RR = RR * RR;
      if (dd <= RR) {
        C.i = i;
        C.j = j;
        if (dd > 0.0)
          C.branch = sqrt(dd);
        else {
          cerr << "dd <= 0 !/n";
          continue;
        }
        C.nx = dx / C.branch;
        C.ny = dy / C.branch;
        C.fn = C.ft = 0.0;
        C.x = 0.5 * (grains[i].x + grains[j].x + C.nx * (Ri - Rj));
        C.y = 0.5 * (grains[i].y + grains[j].y + C.ny * (Ri - Rj));
        neighbors.push_back(C);
      }
    }
  }

  // unpack
  size_t k, kpacked = 0;
  for (k = 0; k < neighbors.size(); ++k) {
    while (kpacked < packed.size() && packed[kpacked].i < neighbors[k].i) ++kpacked;
    if (kpacked == packed.size()) break;

    while (kpacked < packed.size() && packed[kpacked].i == neighbors[k].i && packed[kpacked].j < neighbors[k].j)
      ++kpacked;
    if (kpacked == packed.size()) break;

    if (packed[kpacked].i == neighbors[k].i && packed[kpacked].j == neighbors[k].j) {
      neighbors[k] = packed[kpacked];
      ++kpacked;
    }
  }
}

void Configuration::read(const char* path, int num, int format) {
  char fname[256];
  sprintf(fname, "%s/CONF%04d", path, num);
  ifstream file(fname);
  if (!file) {
    cerr << "A problem occurred when opening the file " << fname << endl;
  }
  grain G;
  int ng = 0;
  grains.clear();
  file >> ng;
  if (format == format_x_y_R_rot) {
    for (int i = 0; i < ng; i++) {
      file >> G.x >> G.y >> G.R >> G.rot;
      grains.push_back(G);
    }
  } else if (format == format_x_y_R_rot_NCC_NCC) {
    for (int i = 0; i < ng; i++) {
      file >> G.x >> G.y >> G.R >> G.rot >> G.zncc1 >> G.zncc2;
      grains.push_back(G);
    }
  } else {
    std::cerr << "@Configuration::read, Unknown format" << std::endl;
    return;
  }
  char line[512];
  char fake;
  file >> fake;  // Maybe the seek position is still at the previous line (file written with fortran)
                 // So, we read a unused (fake) character before ignoring the full line
  file.getline(line, 512);
  file >> Box.x[0] >> Box.y[0];
  file >> Box.x[1] >> Box.y[1];
  file >> Box.x[2] >> Box.y[2];
  file >> Box.x[3] >> Box.y[3];

  file >> fake;
  file.getline(line, 512);

  file >> Strain.Gxx;
  file >> Strain.Gyy;
  file >> Strain.Gxy;
  file >> Strain.Gyx;
  file >> Strain.Exx;
  file >> Strain.Eyy;
  file >> Strain.Exy;
  file >> Strain.Erot;
  file >> Strain.Evol;

  file >> fake;
  file.getline(line, 512);
  file >> scaleFactor;
}

void Configuration::computeDisplacements(Configuration& RefConf) {
  for (size_t i = 0; i < grains.size(); ++i) {
    double xc = grains[i].x;
    double yc = grains[i].y;
    double rot = grains[i].rot;
    double xcRef = RefConf.grains[i].x;
    double ycRef = RefConf.grains[i].y;
    double rotRef = RefConf.grains[i].rot;
    grains[i].dx = (xc - xcRef);
    grains[i].dy = (yc - ycRef);
    grains[i].rot = (rot - rotRef);
  }
}

void Configuration::computeFluctuations(Configuration& RefConf) {

  // Numbering for 1g2e frame:
  //   2 ---- 1
  //   |      |
  //   3 ---- 0
  //
  // Numbering for Bilinear functions:
  //   2 ---- 3
  //   |      |
  //   0 ---- 1

  double x0 = Box.x[3];
  double y0 = Box.y[3];
  double x1 = Box.x[0];
  double y1 = Box.y[0];
  double x2 = Box.x[2];
  double y2 = Box.y[2];
  double x3 = Box.x[1];
  double y3 = Box.y[1];

  double dx0 = x0 - RefConf.Box.x[3];
  double dy0 = y0 - RefConf.Box.y[3];
  double dx1 = x1 - RefConf.Box.x[0];
  double dy1 = y1 - RefConf.Box.y[0];
  double dx2 = x2 - RefConf.Box.x[2];
  double dy2 = y2 - RefConf.Box.y[2];
  double dx3 = x3 - RefConf.Box.x[1];
  double dy3 = y3 - RefConf.Box.y[1];

  double s, t, s2, t2;
  double dxContinuous, dyContinuous;
  for (size_t i = 0; i < grains.size(); ++i) {
    double xc = grains[i].x;
    double yc = grains[i].y;
    double rot = grains[i].rot;
    double xcRef = RefConf.grains[i].x;
    double ycRef = RefConf.grains[i].y;
    double rotRef = RefConf.grains[i].rot;

    grains[i].dx = (xc - xcRef);
    grains[i].dy = (yc - ycRef);
    grains[i].rot = (rot - rotRef);

    inverseBilinear(x0, y0, x1, y1, x2, y2, x3, y3, xc, yc, &s, &t, &s2, &t2);
    bilinear(dx0, dy0, dx1, dy1, dx2, dy2, dx3, dy3, s, t, &(dxContinuous), &(dyContinuous));
    
    grains[i].ddx = (grains[i].dx - dxContinuous);
    grains[i].ddy = (grains[i].dy - dyContinuous);
    grains[i].drot = (grains[i].rot - 0.0);  
  }
}

void Configuration::computeLocalMeanStrains(Configuration& RefConf) {}

double Configuration::getRmin() {
  double Rmin = grains[0].R;
  for (size_t i = 0; i < grains.size(); i++) {
    if (Rmin > grains[i].R) Rmin = grains[i].R;
  }
  return Rmin;
}

double Configuration::getRmax() {
  double Rmax = grains[0].R;
  for (size_t i = 0; i < grains.size(); i++) {
    if (Rmax < grains[i].R) Rmax = grains[i].R;
  }
  return Rmax;
}

double Configuration::getRawOverlay(size_t ineighbor) {
  size_t i = neighbors[ineighbor].i;
  size_t j = neighbors[ineighbor].j;
  double RawOverlay = neighbors[ineighbor].branch - (grains[i].R + grains[j].R);
  return RawOverlay;
}

std::vector<std::pair<size_t, size_t> > buildNetworkCorrespondance(std::vector<neighbor>& A, std::vector<neighbor>& B) {
  std::vector<std::pair<size_t, size_t> > Coupled;

  size_t kA, kB = 0;
  for (kA = 0; kA < A.size(); ++kA) {
    while (kB < B.size() && B[kB].i < A[kA].i) ++kB;
    if (kB == B.size()) break;

    while (kB < B.size() && B[kB].i == A[kA].i && B[kB].j < A[kA].j) ++kB;
    if (kB == B.size()) break;

    if (B[kB].i == A[kA].i && B[kB].j == A[kA].j) {
      // A[kA].ft = B[kB].ft;
      Coupled.push_back(std::pair<size_t, size_t>(kA, kB));
      ++kB;
    }
  }

  return Coupled;
}