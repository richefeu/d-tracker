
struct image_interpolator {
  // In this pur virtual function, x and y are subpixel coordinates
  virtual double getValue(int, double, double) = 0;
};

struct image_interpLinear : public image_interpolator {
  double LinearInterpolate(double p[2], double x) { return ((p[1] - p[0]) * x + p[0]); }

  double getValue(int i, double x, double y) {
    int ix = (int)(floor(x));
    int iy = (int)(floor(y));
    double x01 = x - (double)ix;
    double y01 = y - (double)iy;

    double p[2][2];
    // Il n'y a aucune securite ici -> plantage si x ou y sort des clous
    p[0][0] = image[i][ix][iy];
    p[1][0] = image[i][ix + 1][iy];
    p[0][1] = image[i][ix][iy + 1];
    p[1][1] = image[i][ix + 1][iy + 1];

    double arr[2];
    arr[0] = LinearInterpolate(p[0], y01);
    arr[1] = LinearInterpolate(p[1], y01);

    return LinearInterpolate(arr, x01);
  }
};

// Adapted from: http://www.paulinternet.nl/?page=bicubic
struct image_interpCubic : public image_interpolator {
  double CubicInterpolate(double p[4], double x) {
    double f0 = p[1];
    double f1 = p[2];
    double fp0 = 0.5 * (p[2] - p[0]);
    double fp1 = 0.5 * (p[3] - p[1]);
    double d = f0;
    double c = fp0;
    double b = -3.0 * f0 + 3.0 * f1 - 2.0 * fp0 - fp1;
    double a = 2.0 * f0 - 2.0 * f1 + fp0 + fp1;
    double x2 = x * x;
    return (a * x2 * x + b * x2 + c * x + d);
  }

  double getValue(int i, double x, double y) {
    int ix = (int)(floor(x));
    int iy = (int)(floor(y));
    double y01 = y - (double)iy;
    double x01 = x - (double)ix;
    double p[4][4];
    // Il n'y a aucune securite ici -> plantage si x ou y sort des clous
    p[0][0] = image[i][ix - 1][iy - 1];
    p[1][0] = image[i][ix][iy - 1];
    p[2][0] = image[i][ix + 1][iy - 1];
    p[3][0] = image[i][ix + 2][iy - 1];
    p[0][1] = image[i][ix - 1][iy];
    p[1][1] = image[i][ix][iy];
    p[2][1] = image[i][ix + 1][iy];
    p[3][1] = image[i][ix + 2][iy];
    p[0][2] = image[i][ix - 1][iy + 1];
    p[1][2] = image[i][ix][iy + 1];
    p[2][2] = image[i][ix + 1][iy + 1];
    p[3][2] = image[i][ix + 2][iy + 1];
    p[0][3] = image[i][ix - 1][iy + 2];
    p[1][3] = image[i][ix][iy + 2];
    p[2][3] = image[i][ix + 1][iy + 2];
    p[3][3] = image[i][ix + 2][iy + 2];

    double arr[4];
    arr[0] = CubicInterpolate(p[0], y01);
    arr[1] = CubicInterpolate(p[1], y01);
    arr[2] = CubicInterpolate(p[2], y01);
    arr[3] = CubicInterpolate(p[3], y01);

    return CubicInterpolate(arr, x01);
  }
};

struct image_interpQuintic : public image_interpolator {
  double QuinticInterpolate(double p[4], double x) {
    // values:
    double f0 = p[1];
    double f1 = p[2];
    // First derivatives:
    double fp0 = 0.5 * (p[2] - p[0]);
    double fp1 = 0.5 * (p[3] - p[1]);
    // Second derivatives:
    double fs0 = p[0] - 2.0 * p[1] + p[2];
    double fs1 = p[1] - 2.0 * p[2] + p[3];

    double f = f0;
    double e = fp0;
    double d = 0.5 * fs0;
    double c = 0.5 * (fs1 - 8.0 * fp1 + 20.0 * f1 - 3.0 * fs0 - 12.0 * fp0 - 20.0 * f0);
    double b = 15.0 * f0 + 8.0 * fp0 + 1.5 * fs0 - 15.0 * f1 + 7.0 * fp1 - fs1;
    double a = f1 - b - c - 0.5 * fs0 - fp0 - f0;

    double x2 = x * x;
    double x4 = x2 * x2;
    return (a * x4 * x + b * x4 + c * x2 * x + d * x2 + e * x + f);
  }

  double getValue(int i, double x, double y) {
    int ix = (int)(floor(x));
    int iy = (int)(floor(y));
    double x01 = x - (double)ix;
    double y01 = y - (double)iy;
    double p[4][4];

    // Il n'y a aucune securite ici -> plantage si x ou y sort des clous
    p[0][0] = image[i][ix - 1][iy - 1];
    p[1][0] = image[i][ix][iy - 1];
    p[2][0] = image[i][ix + 1][iy - 1];
    p[3][0] = image[i][ix + 2][iy - 1];
    p[0][1] = image[i][ix - 1][iy];
    p[1][1] = image[i][ix][iy];
    p[2][1] = image[i][ix + 1][iy];
    p[3][1] = image[i][ix + 2][iy];
    p[0][2] = image[i][ix - 1][iy + 1];
    p[1][2] = image[i][ix][iy + 1];
    p[2][2] = image[i][ix + 1][iy + 1];
    p[3][2] = image[i][ix + 2][iy + 1];
    p[0][3] = image[i][ix - 1][iy + 2];
    p[1][3] = image[i][ix][iy + 2];
    p[2][3] = image[i][ix + 1][iy + 2];
    p[3][3] = image[i][ix + 2][iy + 2];

    double arr[4];
    arr[0] = QuinticInterpolate(p[0], y01);
    arr[1] = QuinticInterpolate(p[1], y01);
    arr[2] = QuinticInterpolate(p[2], y01);
    arr[3] = QuinticInterpolate(p[3], y01);

    return QuinticInterpolate(arr, x01);
  }
};

image_interpLinear IMAGE_INTERPOLATOR_LINEAR;
image_interpCubic IMAGE_INTERPOLATOR_CUBIC;
image_interpQuintic IMAGE_INTERPOLATOR_QUINTIC;

image_interpolator *IMAGE_INTERPOLATOR = &IMAGE_INTERPOLATOR_CUBIC;

