/************************************************************************************************/
/*                 DETERMINATION OF PARTICLE CENTERS WITH SUBPIXEL ACCURACY                     */
/************************************************************************************************/

double center_to_minimize_Circ(vector<double> &X) {
  double lx, ly, sum_dist = 0;

  for (size_t i = 0; i < xbound.size(); i++) {
    // undistor(&disto_parameters[0], xd, yd, xu, yu); // TODO...
    lx = xbound[i] - X[0];
    ly = ybound[i] - X[1];
    sum_dist += fabs(sqrt(lx * lx + ly * ly) - X[2]);
  }

  return sum_dist;
}

void find_subpixel_centers() {
  cout << "FIND SUBPIXEL CENTERS" << endl;
  read_image(0, iref, true); // Read image

  // 3 profiles to identify the threshold value of gray
  ofstream prof("profile.txt");
  int yline = (int)(dimy / 4);
  for (int i = 0; i < dimx; i++) {
    prof << i << " " << image[0][i][yline] << endl;
  }
  prof << endl;
  yline = (int)(dimy / 2);
  for (int i = 0; i < dimx; i++) {
    prof << i << " " << image[0][i][yline] << endl;
  }
  prof << endl;
  yline = (int)(3 * dimy / 4);
  for (int i = 0; i < dimx; i++) {
    prof << i << " " << image[0][i][yline] << endl;
  }

  double x, y;
  double x0, y0;
  double x1, y1;
  double x2, y2;
  double x3, y3;
  double gray0, gray1, gray2, gray3;
  double gray, previous_gray;

  ofstream file("recentered.data");
  ofstream logFile("recentered.log.data");
  file << num_grains << endl << endl;
  logFile << num_grains << endl << endl;
  for (int igrain = 0; igrain < num_grains; igrain++) {
    int px0 = grain[igrain].refcoord_xpix;
    int py0 = grain[igrain].refcoord_ypix;
    double r0 = grain[igrain].radius_pix;

    center_parameters[0] = 0.0; // x offset
    center_parameters[1] = 0.0; // y offset
    center_parameters[2] = r0;  // radius

    center_parameters_perturb[0] = subpix_center_xstep; // x offset
    center_parameters_perturb[1] = subpix_center_ystep; // y offset
    center_parameters_perturb[2] = subpix_center_rstep; // radius

    double nx, ny;
    double dr = 0.1;
    xbound.clear();
    ybound.clear();
    for (double delta = 0; delta < 2 * M_PI; delta += M_PI / 180.) {
      nx = cos(delta);
      ny = sin(delta);
      previous_gray = 0;
      for (double r = r0 - subpix_center_dr0; r < r0 + subpix_center_dr0; r += dr) {
        x = px0 + r * nx;
        y = py0 + r * ny;
        x0 = x3 = floor(x);
        y0 = y1 = floor(y);
        x1 = x2 = floor(x + 1.0);
        y2 = y3 = floor(y + 1.0);
        gray0 = image[0][(int)x0][(int)y0];
        gray1 = image[0][(int)x1][(int)y1];
        gray2 = image[0][(int)x2][(int)y2];
        gray3 = image[0][(int)x3][(int)y3];
        // http://en.wikipedia.org/wiki/Bilinear_interpolation
        gray = gray2 * fabs(x - x0) * fabs(y - y0) + gray3 * fabs(x1 - x) * fabs(y - y1) +
               gray0 * fabs(x2 - x) * fabs(y2 - y) + gray1 * fabs(x - x3) * fabs(y3 - y);
        if (gray < subpix_center_threshold && previous_gray > subpix_center_threshold) {
          double a = (gray - previous_gray) / dr;
          double b = gray - a * r;
          double rfit;
          if (fabs(a) < 1e-6)
            rfit = r - dr / 2;
          else
            rfit = (subpix_center_threshold - b) / a;
          xbound.push_back(rfit * nx);
          ybound.push_back(rfit * ny);
          break;
        }
        previous_gray = gray;
      }
    }

    Powell<double(vector<double> &)> powell(center_to_minimize_Circ, 1e-8);
    center_parameters = powell.minimize(center_parameters, center_parameters_perturb);

    cout << "grain number " << igrain << endl;
    cout << "x offset = " << center_parameters[0] << endl;
    cout << "y offset = " << center_parameters[1] << endl;
    cout << "radius   = " << center_parameters[2] << endl;
    cout << endl;

    logFile << px0 << " " << py0 << " " << r0 << " " << center_parameters[0] << " " << center_parameters[1] << " "
            << center_parameters[2] << endl;
    file << px0 + center_parameters[0] << " " << py0 + center_parameters[1] << " " << grain[igrain].refrot << " "
         << center_parameters[2] << endl;
  } // for loop igrain
}
