/************************************************************************************************/
/*                                SYNTHETIC IMAGE GENERATION                                    */
/************************************************************************************************/

// ref for Perlin noise: http://www.dreamincode.net/forums/topic/66480-perlin-noise/

inline double findnoise2(double x, double y) {
  int n = (int)x + (int)y * 57;
  n = (n << 13) ^ n;
  int nn = (n * (n * n * 60493 + 19990303) + 1376312589) & 0x7fffffff;
  return 1.0 - ((double)nn / 1073741824.0);
}

// Cosine interpolation
inline double interpolate(double a, double b, double x) {
  double ft = x * 3.1415927;
  double f = (1.0 - cos(ft)) * 0.5;
  return a * (1.0 - f) + b * f;
}

double noise(double x, double y) {
  double floorx = (double)((int)x); // This is kinda a cheap way to floor a double integer.
  double floory = (double)((int)y);
  double s, t, u, v; // Integer declaration
  s = findnoise2(floorx, floory);
  t = findnoise2(floorx + 1, floory);
  u = findnoise2(floorx, floory + 1); // Get the surrounding pixels to calculate the transition.
  v = findnoise2(floorx + 1, floory + 1);
  double int1 = interpolate(s, t, x - floorx); // Interpolate between the values.
  double int2 = interpolate(u, v, x - floorx); // Here we use x-floorx, to get 1st dimension. Don't mind the x-floorx
                                               // thingie, it's part of the cosine formula.
  return interpolate(int1, int2, y - floory);  // Here we use y-floory, to get the 2nd dimension.
}

// Transformation: first deformation, then rotation, and finally translation
// w and h speak for themselves, zoom will zoom in and out on it (eg zoom = 75).
// P stands for persistence, this controls the roughness of the picture (eg P = 0.5)
void generate_synthetic_images(int w, int h, double zoom, int octaves, double p, double defx, double defy,
                               double transx_pix, double transy_pix, double rot_deg) {
  static int num_image_synth = 0;
  thumbnail synth_im(w, h);

  double rot_rad = rot_deg * (M_PI / 180.);
  double c = cos(rot_rad);
  double s = sin(rot_rad);
  double xc = (double)w * 0.5;
  double yc = (double)h * 0.5;

  double dx, dy;

  // Loop trough all the pixels
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {

      double getnoise = 0;
      for (int a = 0; a < octaves - 1; a++) {   // Loop trough the octaves.
        double frequency = pow(2.0, (double)a); // Increases the frequency
        double amplitude = pow(p, (double)a);   // Decreases the amplitude

        dx = (double)x - xc;
        dy = (double)y - yc;

        double defdx = dx * 1.0 / (1.0 + defx);
        double defdy = dy * 1.0 / (1.0 + defy);
        double decalSecur = 5000.0; // noise n'aime pas les valeurs négatives
        getnoise += noise((decalSecur - transx_pix + xc + (c * defdx - s * defdy)) * (frequency / zoom),
                          (decalSecur - transy_pix + yc + (s * defdx + c * defdy)) * (frequency / zoom)) *
                    amplitude;
      }

      if (getnoise < -1.0)
        getnoise = -1.0;
      else if (getnoise > 1.0)
        getnoise = 1.0;

      synth_im.imThumb[x][y] = UINT16_MAX * (getnoise * 0.5 + 0.5); // Between 0 and 1
    }
  }

  char fname[256];
  sprintf(fname, "synth_oct%d_p%.1f_z%.1f_%d.tif", octaves, (float)p, (float)zoom, num_image_synth++);
  cout << "def = (" << defx << ", " << defy << "),"
       << " trans = (" << transy_pix << ", " << transy_pix << "), rot = " << rot_deg << " deg. -> " << fname
       << endl;
  synth_im.writeTiff(fname, 0.0);
}
