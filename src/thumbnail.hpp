#pragma once

struct RGB {
  unsigned char r, g, b;
  RGB() : r(255), g(255), b(255) {}
  RGB(unsigned char R, unsigned char G, unsigned char B) : r(R), g(G), b(B) {}
  RGB(const RGB &col) : r(col.r), g(col.g), b(col.b) {}
  RGB &operator=(const RGB &other) {
    r = other.r;
    g = other.g;
    b = other.b;
    return *this;
  }
  friend RGB operator+(const RGB &a, const RGB &b) { return RGB(a.r + b.r, a.g + b.g, a.b + b.b); }
  friend RGB operator-(const RGB &a, const RGB &b) { return RGB(a.r - b.r, a.g - b.g, a.b - b.b); }
  friend RGB operator-(const RGB &a) { return RGB(-a.r, -a.g, -a.b); }
  friend RGB operator*(const RGB &a, double k) { return RGB(a.r * k, a.g * k, a.b * k); }
  friend RGB operator*(double k, const RGB &a) { return RGB(a.r * k, a.g * k, a.b * k); }
  friend RGB operator/(const RGB &a, double k) { return RGB(a.r / k, a.g / k, a.b / k); }
  friend void swap(RGB &a, RGB &b) {
    RGB tmp = a;
    a = b;
    b = tmp;
  }
};

class thumbnail {
public:
  vector<vector<uint16_t>> imThumb;
  vector<vector<RGB>> colorThumb; // superimposed colors
  int imW;
  int imH;
  int image_div;

  // Ctor
  thumbnail(vector<vector<vector<uint16_t>>> &image, int index, int imageDiv) : image_div(imageDiv) {
    imW = dimx / image_div; // integer division
    imH = dimy / image_div;
    cout << "Image size: " << dimx << "x" << dimy << endl;
    cout << "Thumbnail size: " << imW << "x" << imH << endl;

    imThumb.resize(imW);
    colorThumb.resize(imW);
    for (int i = 0; i < imW; i++) {
      imThumb[i].resize(imH);
      colorThumb[i].resize(imH);
    }

    // Resize
    int x0 = 0, y0 = 0;
    uint16_t mean_grey;
    double fact;
    int ix = 0, iy = 0;

    uint16_t max_grey = 0;
    uint16_t min_grey = UINT16_MAX;
    for (int xm = image_div; xm < dimx; xm += image_div) {
      y0 = 0;
      iy = 0;
      for (int ym = image_div; ym < dimy; ym += image_div) {
        mean_grey = 0;
        fact = 1.0 / ((xm - x0) * (ym - y0));
        for (int x = x0; x < xm; x++) {
          for (int y = y0; y < ym; y++) {
            mean_grey += nearest(fact * image[index][x][y]);
          }
        }
        imThumb[ix][iy] = mean_grey;

        if (mean_grey < min_grey)
          min_grey = mean_grey;
        if (mean_grey > max_grey)
          max_grey = mean_grey;

        iy++;
        if (iy >= imH)
          iy = imH - 1;

        y0 += image_div;
      }

      x0 += image_div;
      ix++;
      if (ix >= imW)
        ix = imW - 1;
    }

    // Rescale the grey-levels over the whole 16-bit range [0 65535]
    double rescaleFactor = 1.0;
    if ((max_grey - min_grey) > 0)
      rescaleFactor = UINT16_MAX / (max_grey - min_grey);
    for (int y = 0; y < imH; y++) {
      for (int x = 0; x < imW; x++) {
        imThumb[x][y] = (uint16_t)nearest((imThumb[x][y] - min_grey) * rescaleFactor);
      }
    }
  }

  // Ctor for empty (white) image
  thumbnail(int dimx, int dimy) : imW(dimx), imH(dimy), image_div(1) {
    imThumb.resize(imW);
    colorThumb.resize(imW);
    for (int i = 0; i < imW; i++) {
      imThumb[i].resize(imH);
      colorThumb[i].resize(imH);
    }

    for (int y = 0; y < imH; y++) {
      for (int x = 0; x < imW; x++) {
        imThumb[x][y] = UINT16_MAX;
        colorThumb[x][y].r = 0;
        colorThumb[x][y].g = 0;
        colorThumb[x][y].b = 0;
      }
    }
  }

  void writePpm(const char *name, double alpha) {
    ofstream file(name, ios::binary);
    file << "P6" << endl;
    file << imW << " " << imH << endl;
    file << "255" << endl;

    double conv = 255.0 / 65535.0;

    unsigned char v, r, g, b;
    for (int y = 0; y < imH; y++) {
      for (int x = 0; x < imW; x++) {
        v = (unsigned char)nearest(conv * (imThumb[x][y]));

        if (colorThumb[x][y].r == 255 && colorThumb[x][y].g == 255 && colorThumb[x][y].b == 255) {
          r = v;
          g = v;
          b = v;
        } else {
          r = (unsigned char)nearest(alpha * colorThumb[x][y].r + (1.0 - alpha) * v);
          g = (unsigned char)nearest(alpha * colorThumb[x][y].g + (1.0 - alpha) * v);
          b = (unsigned char)nearest(alpha * colorThumb[x][y].b + (1.0 - alpha) * v);
        }

        file.write((const char *)&r, sizeof(unsigned char));
        file.write((const char *)&g, sizeof(unsigned char));
        file.write((const char *)&b, sizeof(unsigned char));
      }
    }
  }

  /// @see https://github.com/abadams/ImageStack/blob/master/src/FileTIFF.cpp
  void writeTiff(const char *name, double alpha) {
    TIFF *output_image;

    // Open the TIFF file
    if ((output_image = TIFFOpen(name, "w")) == NULL) {
      cerr << "Unable to write tif file: " << name << endl;
    }

    // We need to set some values for basic tags before we can add any data
    TIFFSetField(output_image, TIFFTAG_SAMPLESPERPIXEL, 3); // rvb
    TIFFSetField(output_image, TIFFTAG_IMAGEWIDTH, imW);
    TIFFSetField(output_image, TIFFTAG_IMAGELENGTH, imH);
    TIFFSetField(output_image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(output_image, TIFFTAG_ROWSPERSTRIP, 1L);
    TIFFSetField(output_image, TIFFTAG_XRESOLUTION, 1.0);
    TIFFSetField(output_image, TIFFTAG_YRESOLUTION, 1.0);
    TIFFSetField(output_image, TIFFTAG_RESOLUTIONUNIT, 1);
    TIFFSetField(output_image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(output_image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(output_image, TIFFTAG_ORIENTATION, (int)ORIENTATION_TOPLEFT);

    TIFFSetField(output_image, TIFFTAG_BITSPERSAMPLE, 8); // 0..255
    TIFFSetField(output_image, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);

    // Now save the image (line by line)
    double conv = 255.0 / 65535.0;
    vector<unsigned char> m_image_data(imW * 3);
    unsigned char v, r, g, b;
    for (int y = 0; y < imH; y++) {
      for (int x = 0; x < imW; x++) {
        v = (unsigned char)nearest(conv * (imThumb[x][y]));

        if (colorThumb[x][y].r == 255 && colorThumb[x][y].g == 255 && colorThumb[x][y].b == 255) {
          r = v;
          g = v;
          b = v;
        } else {
          r = (unsigned char)nearest(alpha * colorThumb[x][y].r + (1.0 - alpha) * v);
          g = (unsigned char)nearest(alpha * colorThumb[x][y].g + (1.0 - alpha) * v);
          b = (unsigned char)nearest(alpha * colorThumb[x][y].b + (1.0 - alpha) * v);
        }

        m_image_data[x * 3 + 0] = r;
        m_image_data[x * 3 + 1] = g;
        m_image_data[x * 3 + 2] = b;
      }
      TIFFWriteScanline(output_image, &m_image_data[0], y, 1);
    }

    // Close the file
    TIFFClose(output_image);
  }

  void draw_line(int x1, int y1, int x2, int y2, RGB col) {
    // see: http://rosettacode.org/wiki/Bitmap/Bresenham's_line_algorithm#C.2B.2B
    // Bresenham's line algorithm
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if (steep) {
      std::swap(x1, y1);
      std::swap(x2, y2);
    }

    if (x1 > x2) {
      std::swap(x1, x2);
      std::swap(y1, y2);
    }

    const float dx = x2 - x1;
    const float dy = fabs(y2 - y1);

    float error = dx / 2.0f;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for (int x = (int)x1; x < maxX; x++) {
      if (steep) {
        colorThumb[y][x] = col;
      } else {
        colorThumb[x][y] = col;
      }

      error -= dy;
      if (error < 0) {
        y += ystep;
        error += dx;
      }
    }
  }

  void draw_circle(int x_centre, int y_centre, int radius, RGB col) {
    int x = 0, y = radius, m = 5 - 4 * radius;
    while (x <= y) { /// @see http://fr.wikipedia.org/wiki/Algorithme_de_tracé_d'arc_de_cercle_de_Bresenham
      colorThumb[x + x_centre][y + y_centre] = col;
      colorThumb[y + x_centre][x + y_centre] = col;
      colorThumb[-x + x_centre][y + y_centre] = col;
      colorThumb[-y + x_centre][x + y_centre] = col;
      colorThumb[x + x_centre][-y + y_centre] = col;
      colorThumb[y + x_centre][-x + y_centre] = col;
      colorThumb[-x + x_centre][-y + y_centre] = col;
      colorThumb[-y + x_centre][-x + y_centre] = col;
      if (m > 0) {
        y = y - 1;
        m = m - 8 * y;
      }
      x = x + 1;
      m = m + 8 * x + 4;
    }
  }

  void draw_disk(int x_centre, int y_centre, int radius, RGB col) {
    int d = 3 - (2 * radius);
    int x = 0;
    int y = radius;

    while (y >= x) {
      for (int i = x_centre - x; i <= x_centre - x + 2 * x + 1; i++)
        colorThumb[i][y_centre - y] = col;
      for (int i = x_centre - x; i <= x_centre - x + 2 * x + 1; i++)
        colorThumb[i][y_centre + y] = col;
      for (int i = x_centre - y; i <= x_centre - y + 2 * y + 1; i++)
        colorThumb[i][y_centre - x] = col;
      for (int i = x_centre - y; i <= x_centre - y + 2 * y + 1; i++)
        colorThumb[i][y_centre + x] = col;

      if (d < 0)
        d = d + (4 * x) + 6;
      else {
        d = d + 4 * (x - y) + 10;
        y--;
      }

      x++;
    }
  }

  void draw_triangle(int X1, int Y1, int X2, int Y2, int X3, int Y3, RGB Col1, RGB Col2, RGB Col3) {
    // numbering
    int x0, y0;
    int x1, y1;
    int x2, y2;
    int x3, y3;
    RGB col1, col2, col3;

    y3 = std::max(Y1, Y2);
    y3 = std::max(y3, Y3);
    if (y3 == Y1) {
      x3 = X1;
      col3 = Col1;
      y2 = std::min(Y3, Y2);
      if (y2 == Y3) {
        x2 = X3;
        col2 = Col3;
        x1 = X2;
        y1 = Y2;
        col1 = Col2;
      } else {
        x2 = X2;
        col2 = Col2;
        x1 = X3;
        y1 = Y3;
        col1 = Col3;
      }
    } else if (y3 == Y2) {
      x3 = X2;
      col3 = Col2;
      y2 = std::min(Y3, Y1);
      if (y2 == Y3) {
        x2 = X3;
        col2 = Col3;
        x1 = X1;
        y1 = Y1;
        col1 = Col1;
      } else {
        x2 = X1;
        col2 = Col1;
        x1 = X3;
        y1 = Y3;
        col1 = Col3;
      }
    } else {
      x3 = X3;
      col3 = Col3;
      y2 = std::min(Y2, Y1);
      if (y2 == Y2) {
        x2 = X2;
        col2 = Col2;
        x1 = X1;
        y1 = Y1;
        col1 = Col1;
      } else {
        x2 = X1;
        col2 = Col1;
        x1 = X2;
        y1 = Y2;
        col1 = Col2;
      }
    }

    if (y3 == y2)
      return;
    x0 = (x3 - x2) * (y1 - y2) / (y3 - y2) + x2;
    y0 = y1;

    double dx = x0 - x2;
    double dy = y0 - y2;
    double lx = x3 - x2;
    double ly = y3 - y2;
    double l;
    if (lx == 0.0 || ly == 0.0)
      l = 0.0;
    else
      l = sqrt((dx * dx + dy * dy) / (lx * lx + ly * ly));
    RGB col0 = (1.0 - l) * col2 + l * col3;

    // first half
    int H = y0 - y2;
    if (H != 0) {
      double invH = 1.0 / (double)H;
      for (int y = 0; y <= H; y++) {
        double l = (double)y * invH;
        int xd = nearest((1.0 - l) * x2 + l * x0);
        int xe = nearest((1.0 - l) * x2 + l * x1);
        RGB cold = (1.0 - l) * col2 + l * col0;
        RGB cole = (1.0 - l) * col2 + l * col1;
        if (xd > xe) {
          std::swap(xd, xe);
          std::swap(cold, cole);
        }
        double invL = 1.0 / (double)(xe - xd);
        for (int x = xd; x <= xe; x++) {
          double l2 = (x - xd) * invL;
          RGB col = (1.0 - l2) * cold + l2 * cole;
          int yy = y2 + y;
          if (x >= 0 && x < dimx && yy >= 0 && yy < dimy)
            colorThumb[x][yy] = col;
        }
      }
    }

    // second half
    H = y3 - y0;
    if (H != 0) {
      double invH = 1.0 / (double)H;
      for (int y = 1; y <= H; y++) {
        double l = (double)y * invH;
        int xd = nearest((1.0 - l) * x3 + l * x0);
        int xe = nearest((1.0 - l) * x3 + l * x1);
        RGB cold = (1.0 - l) * col3 + l * col0;
        RGB cole = (1.0 - l) * col3 + l * col1;
        if (xd > xe) {
          std::swap(xd, xe);
          std::swap(cold, cole);
        }
        double invL = 1.0 / (double)(xe - xd);
        for (int x = xd; x <= xe; x++) {
          double l2 = (x - xd) * invL;
          RGB col = (1.0 - l2) * cold + l2 * cole;
          int yy = y3 - y;
          if (x >= 0 && x < dimx && yy >= 0 && yy < dimy)
            colorThumb[x][yy] = col;
        }
      }
    }
  }
};