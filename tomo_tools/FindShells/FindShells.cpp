// g++-8 -O3 -Wall -W FindShells.cpp -o FindShells -I ~/Documents/devel/mbox/common
// 
// ======================================

// standard
#include <cstdlib>
#include <tiffio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

// in common
#include "powell_nr3.hpp"

// local
#include "rawFiles.hpp"
#include "Data.hpp"

using namespace std;

bool insideZone(voxel & vox, cylBox & cyl) {
  if (vox.z < cyl.z0) return false;
  if (vox.z > cyl.z1) return false;
  double rate = (double)(vox.z - cyl.z0) / (double)(cyl.z1 - cyl.z0);
  int Rcyl = nearest(cyl.R0 + (cyl.R1 - cyl.R0) * rate); 
  int xcyl =  nearest(cyl.x0 + (cyl.x1 - cyl.x0) * rate);
  int ycyl =  nearest(cyl.y0 + (cyl.y1 - cyl.y0) * rate);
  
  int dx = vox.x - xcyl;
  int dy = vox.y - ycyl;
  int R = nearest(sqrt((double)dx * dx + (double)dy * dy));
  if (R < Rcyl) return true;
  return false;
}

size_t nbVoxSolidinBlock(Data & data, int x, int y , int z) {
  size_t nx = data.dimx;
  size_t nxny = nx * data.dimy;
  int xmin = x - data.halfGridStep;
  int ymin = y - data.halfGridStep;
  int zmin = z - data.halfGridStep;
  int xmax = x + data.halfGridStep;
  int ymax = y + data.halfGridStep;
  int zmax = z + data.halfGridStep;
  size_t nb = 0;
  for (int x = xmin ; x <= xmax ; x++) {
    for (int y = ymin ; y <= ymax ; y++) {
      for (int z = zmin ; z <= zmax ; z++) {
        size_t idx = (size_t)z * nxny + (size_t)y * nx + (size_t)x;
        nb += (size_t) (data.image[idx]);
      }
    } 
  }
  return nb;
}

void maskShell(Data & data, Shell & s) {
  int xmin, ymin, zmin;
  int xmax, ymax, zmax;
  double r2 = s.r_out * s.r_out;
  double h2_4 = 0.25 * s.h * s.h;
  int d = (int)ceil(sqrt(r2 + h2_4)) + 1;
  xmin = s.x - d;
  if (xmin < 0) xmin = 0;
  ymin = s.y - d;
  if (ymin < 0) ymin = 0;
  zmin = s.z - d;
  if (zmin < 0) zmin = 0;
  xmax = s.x + d;
  if (xmax >= (int)data.dimx) xmax = data.dimx - 1;
  ymax = s.y + d;
  if (ymax >= (int)data.dimy) ymax = data.dimy - 1;
  zmax = s.z + d;
  if (zmax >= (int)data.dimz) zmax = data.dimz - 1;
  
  size_t nx = data.dimx;
  size_t nxny = nx * data.dimy;
  double ax = 1.0, ay = 0.0, az = 0.0;
  rotyz(s.roty, s.rotz, ax, ay, az);
  for (int x = xmin ; x <= xmax ; x++) {
    for (int y = ymin ; y <= ymax ; y++) {
      for (int z = zmin ; z <= zmax ; z++) {
        double vx = x - s.x;
        double vy = y - s.y;
        double vz = z - s.z;
        double va = vx * ax + vy * ay + vz * az;
        double tax = vx - va * ax;
        double tay = vy - va * ay;
        double taz = vz - va * az; 
        double ta2 = tax * tax + tay * tay + taz * taz;
        if ( (va * va) < h2_4 && ta2 < r2 ) {
          size_t idx = (size_t)z * nxny + (size_t)y * nx + (size_t)x;
          data.image[idx] = 0;
        } 
      }
    } 
  }
}

struct ErrorFunc {
  Data * data;
  int nAngle;
  int nH; 
  
  std::vector<voxel> voxels;
  
  ErrorFunc(Data *d, int na, int nh): data(d), nAngle(na), nH(nh) { }
  
  void placeVoxels(vector<double> & X) {
    double dAngle = 2.0 * M_PI / (double)nAngle;
    double dh = data->shells_h / (double)nH;
    voxels.clear();
    voxel Vox;
    
    for (double angle = 0.0 ; angle < 2.0 * M_PI ; angle += dAngle) {
      double c = cos(angle);
      double s = sin(angle);
      double yin = data->shells_r_in * c;
      double zin = data->shells_r_in * s;
      double yout = data->shells_r_out * c;
      double zout = data->shells_r_out * s;
      for (double x = -0.5 * data->shells_h ; x < 0.5 * data->shells_h + 0.5 * dh ; x += dh) {
        double vxout = x;
        double vyout = yout;
        double vzout = zout;
        rotyz(X[3], X[4], vxout, vyout, vzout);
        vxout += X[0];
        vyout += X[1];
        vzout += X[2];   
        Vox.x = nearest(vxout);
        Vox.y = nearest(vyout);
        Vox.z = nearest(vzout);
        if (insideZone(Vox, data->zone))
          voxels.push_back(Vox);
        
        double vxin = x;
        double vyin = yin;
        double vzin = zin;
        rotyz(X[3], X[4], vxin, vyin, vzin);
        vxin += X[0];
        vyin += X[1];
        vzin += X[2];   
        Vox.x = nearest(vxin);
        Vox.y = nearest(vyin);
        Vox.z = nearest(vzin);
        if (insideZone(Vox, data->zone))
          voxels.push_back(Vox);
      }
    }   
  }
  
  double operator() (vector<double> & X) {
    placeVoxels(X);
    if (voxels.empty()) return 1.0;
    double rateInside = voxels.size() / ((double)nAngle * (nH + 1));
    if (rateInside < data->insideRate) return 1.0;
    
    double tot = 0.0;
    
    size_t nx = data->dimx;
    size_t nxny = nx * data->dimy;
    size_t nxnynz = nxny * data->dimz;
    for (size_t i = 0 ; i < voxels.size() ; i++) {
      size_t idx = (size_t)voxels[i].z * nxny + (size_t)voxels[i].y * nx + (size_t)voxels[i].x;
      if (idx >= nxnynz) { // This should never happen
        cout << idx << " > " << nxnynz << " !!!!!!!!!!!!!!!!" << endl;
        cout << "x = " << voxels[i].x << endl;
        cout << "y = " << voxels[i].y << endl;
        cout << "z = " << voxels[i].z << endl;
      }
      double val = (double)data->image[idx];
      tot += val;
    }
    return (1.0 - tot / (double)voxels.size());
  }
};


int main(int argc, char const *argv[]) {
	Data data;
		
	if (argc == 2) {
		read_data(argv[1], data);
	}
	else {
		cerr << "Usages: \n" << argv[0] << " InputDataFile" << endl;
		exit(0);
	}
	
  char fname[256];
  sprintf(fname,"%s.raw", data.filename.c_str());
  cout << "read file '" << fname << "'\n";
  data.image = readCompressedRaw(fname, data.dimx, data.dimy, data.dimz);
  cout << "image size: " << data.dimx << " x " << data.dimy << " x " << data.dimz << endl;
  
  vector<double> X(5); 
  vector<double> dX(5);
  dX[0] = data.delta;
  dX[1] = data.delta;
  dX[2] = data.delta;
  dX[3] = 0.01;
  dX[4] = 0.01; 
  
  int step = nearest(2.0 * (double)data.halfGridStep / (double)(data.blockDivNum)); // step inside a block
  
  // generate target blocks
  std::vector<voxel> targets;
  for (int ix = data.halfGridStep ; ix < (int)data.dimx ; ix += 2 * data.halfGridStep) {
    for (int iy = data.halfGridStep ; iy < (int)data.dimy ; iy += 2 * data.halfGridStep) {
      for (int iz = data.halfGridStep ; iz < (int)data.dimy ; iz += 2 * data.halfGridStep) {
        voxel Vox(ix, iy, iz);
        if (insideZone(Vox, data.zone))
          targets.push_back(Vox);
      }
    }
  }
  
  ErrorFunc func(&data, data.nbAngle, data.nbH);
	Powell<ErrorFunc> powell(func);
  size_t nbTargets = targets.size();
  size_t nbVoxMin = (size_t)nearest(3.1415 * (data.shells_r_out * data.shells_r_out 
    - data.shells_r_in * data.shells_r_in) * data.shells_h * 1.2 /*rateVoxMin*/); 
  if (nbTargets > data.MaxNbTargets) nbTargets = data.MaxNbTargets;
  for (size_t itarget = 0 ; itarget < nbTargets; itarget++) {    
    X[3] = X[4] = 0.0;
    
    int xtarget = targets[itarget].x;
    int ytarget = targets[itarget].y;
    int ztarget = targets[itarget].z;
  
    vector<voxel> testPos;
    for (int x = xtarget - data.halfGridStep ; x <= xtarget + data.halfGridStep ; x += step) {
      for (int y = ytarget - data.halfGridStep ; y <= ytarget + data.halfGridStep ; y += step) {
        for (int z = ztarget - data.halfGridStep ; z <= ztarget + data.halfGridStep ; z += step) {
          if ( x >= 0 && x < (int)data.dimx 
            && y >= 0 && y < (int)data.dimy 
            && z >= 0 && z < (int)data.dimz 
          ) testPos.push_back(voxel(x, y, z));
        }
      }
    }
  
    int found = 0;
    size_t best_trial = 0;
    double best_match = 0.0;
    for (size_t trial = 0 ; trial < testPos.size() ; trial++) {
      X[0] = (double)testPos[trial].x;
      X[1] = (double)testPos[trial].y;
      X[2] = (double)testPos[trial].z;
      X[3] = X[4] = 0.0;
        
      X = powell.minimize(X, dX);
      double match = 1.0 - powell.fret;
      if (match > best_match) { 
        best_match = match;
        best_trial = trial;
      }
    	std::cout << "Target " << itarget + 1 << "/" << nbTargets  
        << ", Trial " << trial + 1 << "/" << testPos.size() << ", match = " << match << std::endl;
    
      if (powell.fret < (1.0 - data.matchScore)) { // fret should actually be strictly nul
        std::cout << "      > One shell found! :)\n";
        found += 1;
        
        Shell s;
        s.x = nearest(X[0]);
        s.y = nearest(X[1]);
        s.z = nearest(X[2]);
        s.roty = X[3];
        s.rotz = X[4];
        s.r_in = data.shells_r_in;
        s.r_out = data.shells_r_out;
        s.h = data.shells_h;
        s.best_match = best_match;
        data.shells.push_back(s);
        maskShell(data, s);
        std::cout << s.x << " " << s.y << " " << s.z << "  "
          << s.roty << " " << s.rotz << "  "
          << s.r_in << " " << s.r_out << " " << s.h << " " << best_match << std::endl;
        
        // Here we test if a given amount of solid still remain in the treated block
        // before leaving
        size_t nbVox = nbVoxSolidinBlock(data, 
          testPos[trial].x, testPos[trial].y, testPos[trial].y);
        if (nbVox <= nbVoxMin) {
          break;
        } 
        std::cout << "continue in the same block!\n";
      }
    }
  
    if (found > 0) {
      std::cout << "      > Target with " << found << " matches :)\n";
    }
    else {
      std::cout << "      > Target without match! :(\n";
      if (best_match > 0.95/*allowedMatch*/) {
        X[0] = (double)testPos[best_trial].x;
        X[1] = (double)testPos[best_trial].y;
        X[2] = (double)testPos[best_trial].z;
        X[3] = X[4] = 0.0;
        
        X = powell.minimize(X, dX);
        std::cout << "      > One shell added with match = " << 1.0 - powell.fret << '\n';
        Shell s;
        s.x = nearest(X[0]);
        s.y = nearest(X[1]);
        s.z = nearest(X[2]);
        s.roty = X[3];
        s.rotz = X[4];
        s.r_in = data.shells_r_in;
        s.r_out = data.shells_r_out;
        s.h = data.shells_h;
        s.best_match = best_match;
        data.shells.push_back(s);
        maskShell(data, s);
        std::cout << s.x << " " << s.y << " " << s.z << "  "
          << s.roty << " " << s.rotz << "  "
          << s.r_in << " " << s.r_out << " " << s.h << " " << s.best_match << std::endl;
      }
    }
        
  } // end for-loop itarget
  
  // Save
  sprintf(fname,"%s.txt", data.filename.c_str());
  std::cout << "Save shell positions in " << fname << std::endl;
  std::ofstream posFile(fname);
  posFile << data.shells.size() << std::endl;
  for (size_t i = 0 ; i < data.shells.size() ; i++) {
    posFile << data.shells[i].x << " " << data.shells[i].y << " " << data.shells[i].z << "  "
        << data.shells[i].roty << " " << data.shells[i].rotz << "  "
        << data.shells[i].r_in << " " << data.shells[i].r_out << " " << data.shells[i].h << " "
        << data.shells[i].best_match << std::endl;
  }
   
  sprintf(fname,"%s_masked.raw", data.filename.c_str());
  std::cout << "Save masked raw-file in " << fname << std::endl;
  saveCompressedRaw(fname, data.image, data.dimx, data.dimy, data.dimz);
    
  return 0;
}



