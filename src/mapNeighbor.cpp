#include "tracker.cpp"
#include <iostream>

int main() {

  read_image(0, "/Users/vrichefeu/Documents/devel/tracker2/examples/ImageData/TEST01.tiff", true);
#if 0
  for(int i = 0 ;i<=100;i++) {
     double xsub = (i/10.) + 147 - 4;
    for (int j = 0 ; j <= 100;j++) {
      double ysub = (j/10.) + 257 - 4;
      double val = IMAGE_INTERPOLATOR_CUBIC.getValue(0, xsub, ysub);
      //std::cout << xsub << " " << ysub << " " << val << "\n";
      if (i%10 == 0 && j%10 == 0) {
        std::cout << xsub << " " << ysub << " " << val << "\n";
      }
    }
    std::cout << "\n";
  }
  return 0;
#endif
  read_image(1, "/Users/vrichefeu/Documents/devel/tracker2/examples/ImageData/TEST02.tiff", false);

  // reference position
  int x0 = 147; // -> 156 (dx = 9)
  int y0 = 257; // -> 267 (dy = 10)
  int x1 = 156;
  int y1 = 267;
  
  //int x0 = 169;
  //int y0 = 321;

  int h = 16; // pattern half size
  
  // convert NCCmap_g1_5x5.tif -crop 5x5+155+264 -interpolate Integer -filter point NCCmap_g1_5x5_zoom.tif
  // convert NCCmap_g1_9x9.tif -crop 9x9+152+263 -interpolate Integer -filter point NCCmap_g1_9x9_zoom.tif
  // convert NCCmap_g1_17x17.tif -crop 17x17+148+259 -interpolate Integer -filter point NCCmap_g1_17x17_zoom.tif
  // convert NCCmap_g1_33x33.tif -crop 33x33+140+251 -interpolate Integer -filter point NCCmap_g1_33x33_zoom.tif
  // convert NCCmap_g1_65x65.tif -crop 65x65+124+235 -interpolate Integer -filter point NCCmap_g1_65x65_zoom.tif

  // window
  int right = dimx - x0 - h - 1;
  int left = x0 - h - 1;
  int up = y0 - h - 1;
  int down = dimy - y0 - h - 1;

  cout << "gauche " << x0 - left << '\n';
  cout << "droite " << x0 + right << '\n';
  cout << "haut   " << x0 + up << '\n';
  cout << "bas    " << x0 - down << '\n';
  cout << right + left << "x" << up + down << "+" << x0 - left << "+" << x0 - down << '\n';

  double nb = (2. * h + 1.) * (2. * h + 1.);

  // precomputation
  double C0C0 = 0.0;
  double mean0 = 0.0;
  double diff = 0.0;
  for (int ix = x0 - h; ix <= x0 + h; ix++) {
    for (int iy = y0 - h; iy <= y0 + h; iy++) {
      mean0 += (double)image[0][ix][iy];
    }
  }
  mean0 /= nb;
  for (int ix = x0 - h; ix <= x0 + h; ix++) {
    for (int iy = y0 - h; iy <= y0 + h; iy++) {
      diff = (double)image[0][ix][iy] - mean0;
      C0C0 += diff * diff;
    }
  }
  
  double NCCMax = -1.;

  // build the map
  double C1C1 = 0.0, C0C1 = 0.0;
  double mean1 = 0.0;
  vector<vector<double>> mapNCC;
  size_t nx = left + right + 1;
  size_t ny = up + down + 1;
  mapNCC.resize(nx);
  for (size_t i = 0; i < nx; i++)
    mapNCC[i].resize(ny);

  size_t xmap = 0, ymap = 0;
  double NCC;
  for (int px = x0 - left; px <= x0 + right; px++) {
    ymap = 0;
    for (int py = y0 - up; py <= y0 + down; py++) {

      C0C1 = C1C1 = mean1 = 0.0;
      for (int ix = px - h; ix <= px + h; ix++) {
        for (int iy = py - h; iy <= py + h; iy++) {
          mean1 += (double)image[1][ix][iy];
        }
      }
      mean1 /= nb;
      for (int ix = -h; ix <= h; ix++) {
        for (int iy = -h; iy <= h; iy++) {
          diff = (double)image[1][px + ix][py + iy] - mean1;
          C0C1 += ((double)image[0][x0 + ix][y0 + iy] - mean0) * diff;
          C1C1 += diff * diff;
        }
      }

      NCC = C0C1 / sqrt(C0C0 * C1C1);
      if (NCC > NCCMax) {
        NCCMax = NCC;
        x1 = xmap;
        y1 = ymap;
      }
      mapNCC[xmap][ymap] = NCC;
      ymap++;
    }
    xmap++;
  }

  // save it as a pgm image
  ColorTable ct;
  ct.setTableID(4); // rainbow
  ct.setMinMax(0.0, 1.0);
  ct.SavePpm("colorGradient.ppm");
  colorRGBA col;

  thumbnail colIm(image, 1, 1);

  for (size_t y = 0; y < ny; y++) {
    for (size_t x = 0; x < nx; x++) {
      ct.getRGB((float)mapNCC[x][y], &col);

      colIm.colorThumb[x0 - left + x][y0 - up + y].r = col.r;
      colIm.colorThumb[x0 - left + x][y0 - up + y].g = col.g;
      colIm.colorThumb[x0 - left + x][y0 - up + y].b = col.b;
      // cout <<  col.r << " " << col.g << " " << col.b << endl;
    }
  }

  colIm.writeTiff("NCCmap.tif", 0.5);
  
  cout << "position = (" << x1 << ", " << y1 << ")\n";
  
  std::ofstream ffx("NCC-lignex_h16.txt");
  for (size_t x = 0; x < nx; x++) {
    ffx << (double)x-(double)x0 << " " << mapNCC[x][y1] << '\n';
  }

  std::ofstream ffy("NCC-ligney_h16.txt");
  for (size_t y = 0; y < ny; y++) {
    ffy << (double)y-(double)y0 << " " << mapNCC[x1][y] << '\n';
  }
  

  return 0;
}