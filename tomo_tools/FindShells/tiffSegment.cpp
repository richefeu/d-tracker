// g++-8 -O3 -Wall tiffSegment.cpp -o tiffSegment -ltiff

#include <cstdlib>
#include <tiffio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

#include "rawFiles.hpp"

using namespace std;

// Read tiff image (16-bits Gray levels)
void readTiff(const char * name, vector< vector<unsigned short> > & image, size_t & dimx, size_t & dimy)
{
	cout << "-> Read " << name << endl;
	TIFF* tif = TIFFOpen(name, "r");
	
	if (!tif) {
		cerr << "Cannot open tiff file named '" << name << "'" << endl;
		exit(0);
	}
	uint32 w, h;
	size_t npixels;
	uint32* raster;
 
	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
	cout << w << "x" << h << endl;
	npixels = w * h;
 
	raster = (uint32*) _TIFFmalloc(npixels * sizeof (uint32));
	if (raster != NULL) {
		if (!TIFFReadRGBAImage(tif, w, h, raster, 0)) {
			cerr << "Cannot read tiff file named '" << name << "'" << endl;
			exit(0);
		}
	}

	dimx = (size_t)w;
	dimy = (size_t)h;

	// Reserve memory for image
	if (!image.empty()) image.clear();
	image.resize(dimx);
		
	for (size_t ix = 0 ; ix < dimx ; ix++) {
		image[ix].resize(dimy);
	}
	
	unsigned short r, g, b;
	size_t dimxy = dimx * dimy;
	for (size_t y = 0 ; y < dimy ; ++y) {
		size_t shift = y * dimx;
		for (size_t x = 0 ; x < dimx ; ++x) {
			r = (unsigned short) raster[dimxy - (shift + (dimx - x))    ];
			g = (unsigned short) raster[dimxy - (shift + (dimx - x)) + 1];
			b = (unsigned short) raster[dimxy - (shift + (dimx - x)) + 2];

			image[x][y] = (unsigned short)(0.299 * r + 0.587 * g + 0.114 * b);
		}
	}
		
	_TIFFfree(raster);
	TIFFClose(tif);
}

unsigned char conv(unsigned short Treshold, unsigned short value) {
  if (value < Treshold) return 0; 
  return 1;
};

// Usage:
// tiffSegment fileNameTemplate imFirst imLast GreyCut
// 
// At the end, the components x and y are the ones of the pictures,
// and the z is the file number
//
// ./tiffSegment ~/Desktop/SlicesY/slice%05d.tif 62 1344
int main (int argc, char const *argv[]) {	
		
	if (argc != 5) {
		cerr << "Usages: \n" << argv[0] << " fileNameTemplate imFirst imLast GreyCut" << endl;
		exit(0);
	}
	
  size_t imFirst = (size_t)atoi(argv[2]);
  size_t imLast = (size_t)atoi(argv[3]);
  short int GreyCut = atoi(argv[4]);
  
  vector<vector<unsigned short> > image;
  size_t dimx, dimy;
  
  char fname[256];
  sprintf(fname, argv[1], imFirst);
  readTiff(fname, image, dimx, dimy);
  size_t dimz = imLast - imFirst + 1;
  cout << "dimx = " << dimx << endl;
  cout << "dimy = " << dimy << endl;
  cout << "dimz = " << dimz << endl;
  size_t ntot = (size_t)dimx * (size_t)dimy * (size_t)dimz;
  cout << "ntot = " << ntot << endl;
  unsigned char* img = new unsigned char[ntot];
  
  size_t id = 0;
  for (size_t i = imFirst ; i <= imLast ; i++) {
    sprintf(fname, argv[1], i);
    readTiff(fname, image, dimx, dimy);
    
    for (size_t y = 0 ; y < dimy ; y++) {
      for (size_t x = 0 ; x < dimx ; x++) {
        img[id++] = conv(GreyCut, image[x][y]);
      }
    }
  }
  cout << "All files have been read" << endl;
  
  cout << "Now, compress and save in a single file named '3Dscan-compr.raw'" << endl;
  saveCompressedRaw("3Dscan-compr.raw", img, dimx, dimy, dimz);
  
  //saveVtrBinary("3Dscan-compr.vtr", img, dimx, dimy, dimz, 1.0);
  cout << "done." << endl;
	return 0;
}
