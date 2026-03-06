#ifndef DATA_HPP_DCE3770B
#define DATA_HPP_DCE3770B

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// Carrot
struct cylBox {
  int x0, y0, z0, R0; // in pixels
  int x1, y1, z1, R1;
}; 

// Not broken shell
struct Shell {
	int x, y, z;
  double roty, rotz; // the non rotated axis is x, and the rotations apply in the order "y then z"
  int r_in, r_out, h; // all integers expresse in pixel units
  double best_match;
  
	Shell(): x(0), y(0), z(0), roty(0.0), rotz(0.0), r_in(0), r_out(0), h(0), best_match(0.0) { }
	Shell(int X, int Y, int Z, int R_IN, int R_OUT, int H): 
    x(X), y(Y), z(Z), roty(0.0), rotz(0.0),
    r_in(R_IN), r_out(R_OUT), h(H), best_match(0.0) { }
  void getAxis(double & nx, double & ny, double & nz) {
    // rotation /y
    double c1 = cos(roty);
    double s1 = sin(roty);
    // rotation /z
    double c2 = cos(rotz);
    double s2 = sin(rotz);
    nx = c2 * c1;
    ny = s2 * c1;
    nz = -s1 ;    
  }
};

struct voxel {
  int x, y, z;
  voxel(): x(0), y(0), z(0) { }
  voxel(int X, int Y, int Z): x(X), y(Y), z(Z) { }
};


void rotyz(double roty, double rotz, double & vx, double & vy, double & vz) {
  // rotation /y
  double c1 = cos(roty);
  double s1 = sin(roty);
  double vvx = c1 * vx + s1 * vz;
  double vvy = vy;
  double vvz = -s1 * vx + c1 * vz;
  // rotation /z
  double c2 = cos(rotz);
  double s2 = sin(rotz);
  vx = c2 * vvx - s2 * vvy;
  vy = s2 * vvx + c2 * vvy;
  vz = vvz;
}

inline int nearest(double x) {
	return (int)floor(x + 0.5);
}

struct Data {
	unsigned char* image; // image (already segmented)
	size_t dimx, dimy, dimz;
  cylBox zone;
	
	std::string filename;
  int shells_r_in, shells_r_out, shells_h;
  std::vector<Shell> shells;
    
  // Algo
  int blockDivNum;
  int nbAngle;
  int nbH;
  int halfGridStep; // half size of sub-samples (blocks)
  int delta;
  size_t MaxNbTargets;
  double insideRate;
  double matchScore;
  
  Data() {
    shells_r_in = 65;
    shells_r_out = 75;
    shells_h = 150;
    nbAngle = 24;
    nbH = 40;
    insideRate = 0.5;
    matchScore = 0.99;
    blockDivNum = 3;
    halfGridStep = 80;
    MaxNbTargets = 100000;
  }
};

void read_data(const char * name, Data & data)
{
	std::ifstream dataFile(name);
	if (!dataFile) {
		std::cerr << "Cannot open file " << name << std::endl;
		exit(0);
	}

	std::string token;
	dataFile >> token;

	while (dataFile) {

		if (token[0] == '!') getline(dataFile, token);
		else if (token == "filename") {
			dataFile >> data.filename;
		}
		else if (token == "blockDivNum") {
			dataFile >> data.blockDivNum;
		}
		else if (token == "nbAngle") {
			dataFile >> data.nbAngle;
		}
		else if (token == "nbH") {
			dataFile >> data.nbH;
		}
		else if (token == "halfGridStep") {
			dataFile >> data.halfGridStep;
		}
		else if (token == "delta") {
			dataFile >> data.delta;
		}
		else if (token == "MaxNbTargets") {
			dataFile >> data.MaxNbTargets;
		}
		else if (token == "insideRate") {
			dataFile >> data.insideRate;
		}
		else if (token == "matchScore") {
			dataFile >> data.matchScore;
		}
		else if (token == "shells.r_out") {
			dataFile >> data.shells_r_out;
		}
		else if (token == "shells.r_in") {
			dataFile >> data.shells_r_in;
		}
		else if (token == "shells.h") {
			dataFile >> data.shells_h;
		}
		else if (token == "zone.x0") {
			dataFile >> data.zone.x0;
		}
		else if (token == "zone.y0") {
			dataFile >> data.zone.y0;
		}
		else if (token == "zone.z0") {
			dataFile >> data.zone.z0;
		}
		else if (token == "zone.x1") {
			dataFile >> data.zone.x1;
		}
		else if (token == "zone.y1") {
			dataFile >> data.zone.y1;
		}
		else if (token == "zone.z1") {
			dataFile >> data.zone.z1;
		}
		else if (token == "zone.R0") {
			dataFile >> data.zone.R0;
		}
		else if (token == "zone.R1") {
			dataFile >> data.zone.R1;
		}
		else {
			fprintf(stdout, "Unknown token: %s\n", token.c_str());
		}

		dataFile >> token;
	}
}


#endif /* end of include guard: DATA_HPP_DCE3770B */
