#ifndef RAWFILES_HPP
#define RAWFILES_HPP

#include <iostream>
#include <fstream>
#include <cmath>

// This can be used to 'compress' the data in a file
struct voxelType {
  // 4 * two-bit field,
  // allowed values are 0...3 for each value
  unsigned char v1:2, v2:2, v3:2, v4:2;
};

struct comprVoxelType {
  unsigned short nb:14, v:2; // nb = 0...16383, v = 0...3
};

/*
struct IMG {
  unsigned char* img;
  size_t nx;
  size_t ny;
  size_t nz;
  size_t ntot;
};

struct box {
  size_t xmin, xmax;
  size_t ymin, ymax;
  size_t zmin, zmax;
};
*/

void saveRaw(const char *name, unsigned char* img, size_t nx, size_t ny, size_t nz) {
  std::ofstream outFile;
  outFile.open(name, std::ios::out | std::ios::binary);
  outFile.write((char*)&nx, sizeof(size_t));
  outFile.write((char*)&ny, sizeof(size_t));
  outFile.write((char*)&nz, sizeof(size_t));
  size_t ntot = nx * ny * nz;
  outFile.write((char*)&img, ntot * sizeof(unsigned char));
  outFile.close();
}


/*
unsigned char* readRaw(const char *name, size_t &nx, size_t &ny, size_t &nz) {
  
}
*/

void savePackedRaw(const char *name, unsigned char* img, size_t nx, size_t ny, size_t nz) {
  std::ofstream outFileCompr;
  outFileCompr.open(name, std::ios::out | std::ios::binary);
  outFileCompr.write((char*)&nx, sizeof(size_t));
  outFileCompr.write((char*)&ny, sizeof(size_t));
  outFileCompr.write((char*)&nz, sizeof(size_t));
  size_t ntot = nx * ny * nz;
  for (size_t i = 0 ; i < ntot ; i += 4) {
    voxelType vox;
    vox.v1 = img[i];
    size_t inext = i + 1;
    vox.v2 = (inext < ntot) ? img[inext] : 0;
    inext++;
    vox.v3 = (inext < ntot) ? img[inext] : 0;
    inext++;
    vox.v4 = (inext < ntot) ? img[inext] : 0;
    outFileCompr.write((char*)&vox, sizeof(voxelType));
  }
  outFileCompr.close();
}

unsigned char* readPackedRaw(const char *name, size_t &nx, size_t &ny, size_t &nz) {
  std::ifstream file(name, std::ios::in | std::ios::binary);
  file.read((char*)&nx, sizeof(size_t));
  file.read((char*)&ny, sizeof(size_t));
  file.read((char*)&nz, sizeof(size_t));
  size_t ntot = nx * ny * nz;
  unsigned char* img = new unsigned char[ntot];
  for (size_t i = 0 ; i < ntot ; i += 4) {
    voxelType vox;
    file.read((char*)&vox, sizeof(voxelType));
    img[i] = (int)vox.v1;
    size_t inext = i + 1;
    if (inext < ntot) img[inext] = (int)vox.v2;
    inext++;
    if (inext < ntot) img[inext] = (int)vox.v3;
    inext++;
    if (inext < ntot) img[inext] = (int)vox.v4;
  }
  file.close();
  return img;
}

void saveCompressedRaw(const char *name, unsigned char* img, size_t nx, size_t ny, size_t nz) {
  std::ofstream outFileCompr;
  outFileCompr.open(name, std::ios::out | std::ios::binary);
  outFileCompr.write((char*)&nx, sizeof(size_t));
  outFileCompr.write((char*)&ny, sizeof(size_t));
  outFileCompr.write((char*)&nz, sizeof(size_t));
  size_t ntot = nx * ny * nz;
  
  size_t nc = 0;
  size_t i0 = 0;
  size_t blockMax = 16383;
  short v;
  size_t nb;
  
  size_t i = 0;
  while (i < ntot) {
    v = img[i0];
    nb = 0;
    while (v == img[i] && i < ntot && nb < blockMax) {
      i++;
      nb = i - i0;
    }
    
    // store
    comprVoxelType blk;
    blk.nb = (short)nb;
    blk.v = v;
    outFileCompr.write((char*)&blk, sizeof(comprVoxelType));
    
    i0 = i;
    nc++; 
  }
  outFileCompr.close();
}

unsigned char* readCompressedRaw(const char *name, size_t &nx, size_t &ny, size_t &nz) {
  std::ifstream file(name, std::ios::in | std::ios::binary);
  file.read((char*)&nx, sizeof(size_t));
  file.read((char*)&ny, sizeof(size_t));
  file.read((char*)&nz, sizeof(size_t));
  size_t ntot = nx * ny * nz;
  //std::cout << "nx = " << nx << std::endl;
  //std::cout << "ny = " << ny << std::endl;
  //std::cout << "nz = " << nz << std::endl;
  //std::cout << "ntot = " << ntot << std::endl;
  unsigned char* img = new unsigned char[ntot];
  size_t i0 = 0;
  while (file.good()) {
    comprVoxelType blk;
    file.read((char*)&blk, sizeof(comprVoxelType));
    size_t i1 = i0 + (size_t)blk.nb;
    for (size_t i = i0 ; i < i1 ; i++) img[i] = (int)blk.v;
    i0 = i1;
    if (i0 == ntot) break;
  }
  file.close();

  return img;
}

void count(unsigned char* img,
           size_t nx, size_t ny, size_t nz,
           size_t x0, size_t y0, size_t z0, 
           size_t x1, size_t y1, size_t z1,
           size_t &Vg, size_t &Vl, size_t &Vs, size_t &Vw) {
  
  size_t ntot = nx * ny * nz;
  size_t nxny = nx * ny;
  if (z0 * nxny + y0 * nx + x0 > ntot) {
    std::cout << "out of range 0\n";
    std::cout << x0 << " " << y0 << " " << z0 << "\n";
  }
  if (z1 * nxny + y1 * nx + x1 > ntot) {
    std::cout << "out of range 1\n";
    std::cout << x1 << " " << y1 << " " << z1 << "\n";
  }
  Vg = Vl = Vs = Vw = 0;
  for (size_t z = z0 ; z < z1 ; z++) {  
    for (size_t y = y0 ; y < y1 ; y++) {
      for (size_t x = x0 ; x < x1 ; x++) {
        size_t idx = z * nxny + y * nx + x;
        if (img[idx] == 0) Vg++;
        else if (img[idx] == 1) Vl++;
        else if (img[idx] == 2) Vs++;
        else if (img[idx] == 3) Vw++;
      }
    }
  }
}


void write_size_t(FILE *fp, size_t val)
{
	fwrite(&val, sizeof(size_t), 1, fp);
}

void write_uchar_t(FILE *fp, unsigned char val)
{
	fwrite(&val, sizeof(unsigned char), 1, fp);
}

void write_float(FILE *fp, float val)
{
	fwrite(&val, sizeof(float), 1, fp);
}


void saveVtrBinary(const char *name, unsigned char* img, 
                   size_t nx, size_t ny, size_t nz, float pasxyz)
{
	FILE * sortie;
	sortie = fopen(name, "wb");

  short int word = 0x0001;
	char *byte = (char*)&word;
	bool big_endian = (byte[0] ? false : true);

	size_t size = 0;
	if (big_endian == true) fprintf(sortie, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"BigEndian\" >\n");
	else fprintf(sortie, "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\" >\n");
  
	fprintf(sortie, "   <RectilinearGrid  WholeExtent=\"0 %zu 0 %zu 0 %zu\">\n", nx - 1, ny - 1, nz - 1);
	fprintf(sortie, "        <Piece Extent=\"0 %zu 0 %zu 0 %zu\">\n", nx - 1, ny - 1, nz - 1);
	fprintf(sortie, "            <Coordinates>\n");
	fprintf(sortie, "                <DataArray type=\"Float32\" Name=\"X\" format=\"appended\" offset=\"%zu\">\n", size);
	size += sizeof(size_t) + nx * sizeof(float);
	fprintf(sortie, "                </DataArray>\n");
	fprintf(sortie, "                <DataArray type=\"Float32\" Name=\"Y\" format=\"appended\" offset=\"%zu\">\n", size);
	size += sizeof(size_t) + ny * sizeof(float);
	fprintf(sortie, "                </DataArray>\n");
	fprintf(sortie, "                <DataArray type=\"Float32\" Name=\"Z\" format=\"appended\" offset=\"%zu\">\n", size);
	size += sizeof(size_t) + nz * sizeof(float);
	fprintf(sortie, "                </DataArray>\n");
	fprintf(sortie, "            </Coordinates>\n");

	fprintf(sortie, "            <PointData Scalars=\"V\" >\n");
	fprintf(sortie, "                <DataArray type=\"UInt8\" Name=\"V\" format=\"appended\" offset=\"%zu\">\n", size);
	size += sizeof(size_t) + sizeof(unsigned char) * nx * ny * nz;
	fprintf(sortie, "                </DataArray>\n");
	fprintf(sortie, "            </PointData>\n");

	fprintf(sortie, "        </Piece>\n");
	fprintf(sortie, "    </RectilinearGrid>\n");
	fprintf(sortie, "    <AppendedData encoding=\"raw\">\n");
	fprintf(sortie, "_");

  // X
	write_size_t(sortie, nx * sizeof(float));
	for (size_t i = 0; i <= nx - 1; i++) {
		write_float(sortie, i * pasxyz);
	}
  
  // Y
	write_size_t(sortie, ny * sizeof(float));
	for (size_t i = 0; i <= ny - 1; i++) {
		write_float(sortie, i * pasxyz);
	}

  // Z
	write_size_t(sortie, nz * sizeof(float));
	for (size_t i = 0; i <= nz - 1; i++) {
		write_float(sortie, i * pasxyz);
	}

	// Values
  size_t nxny = nx * ny;
	write_size_t(sortie, nx * ny * nz * sizeof(unsigned char));
	for (size_t z = 0; z < nz; z++) {
		for (size_t y = 0; y < ny; y++) {
			for (size_t x = 0; x < nx; x++) {
				write_uchar_t(sortie, (unsigned char) img[z * nxny + y * nx + x]);
			}
		}
	}

	fprintf(sortie, "    </AppendedData>\n");
	fprintf(sortie, "</VTKFile>\n");
	fclose(sortie);
}


#endif /* end of include guard: RAWFILES_HPP */
