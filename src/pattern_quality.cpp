/************************************************************************************************/
/*                     DETERMINATION OF PATTERN QUALITY FOR CORRELATIONS                        */
/************************************************************************************************/

double autocorrelation(int x0, int y0, double x1, double y1, int half) {
  double g0, g1;
  double I0mean = 0.0, I1mean = 0.0;
  double sum01 = 0.0, sum00 = 0.0, sum11 = 0.0;

  for (int i = -half; i <= half; i++) {
    for (int j = -half; j <= half; j++) {
      I0mean += image[0][x0 + i][y0 + j];
      I1mean += IMAGE_INTERPOLATOR->getValue(0, x1 + (double)i, y1 + (double)j);
    }
  }
  double nb = (double)(2 * half + 1);
  nb = nb * nb;
  I0mean /= nb;
  I1mean /= nb;

  for (int i = -half; i <= half; i++) {
    for (int j = -half; j <= half; j++) {
      g0 = image[0][x0 + i][y0 + j] - I0mean;
      g1 = IMAGE_INTERPOLATOR->getValue(0, x1 + (double)i, y1 + (double)j) - I1mean;
      sum01 += (g0 * g1);
      sum00 += (g0 * g0);
      sum11 += (g1 * g1);
    }
  }

  return (sum01 / sqrt(sum00 * sum11));
}

void pattern_quality() {
  cout << "PATTERN QUALITY" << endl;
  read_image(0, iref, true); // Read image
  const int nbins = 100;
  int bin = 0;
  ofstream file("autocor.txt");
  vector<vector<double>> autoCor(num_grains);
  for (int i = 0; i < num_grains; i++)
    autoCor[i].resize(nbins);
  vector<vector<double>> population(num_grains);
  for (int i = 0; i < num_grains; i++)
    population[i].resize(nbins);
  double step = 0.05;
  double dist_max = sqrt(2.0 * 10.0 * 10.0);

  for (int igrain = 0; igrain < num_grains; igrain++) {

    cout << "Point " << igrain << endl;

    int x0 = grain[igrain].refcoord_xpix, y0 = grain[igrain].refcoord_ypix; // origin for auto-correlations
    double xmin = x0 - 10.0, xmax = x0 + 10.0;
    double ymin = y0 - 10.0, ymax = y0 + 10.0;

    for (int i = 0; i < nbins; i++) {
      autoCor[igrain][i] = 0.0;
      population[igrain][i] = 0;
    }

    vector<double> x1;
    for (double X1 = xmin; X1 <= xmax; X1 += step)
      x1.push_back(X1);

#pragma omp parallel for private(bin)
    for (size_t ix1 = 0; ix1 < x1.size(); ++ix1) {
      for (double y1 = ymin; y1 <= ymax; y1 += step) {
        double C = autocorrelation(x0, y0, x1[ix1], y1, halfPatternQual);
        double dist = sqrt((x1[ix1] - x0) * (x1[ix1] - x0) + (y1 - y0) * (y1 - y0));
        bin = (int)floor(nbins * dist / dist_max);
        if (bin < nbins) {
          autoCor[igrain][bin] += C;
          population[igrain][bin] += 1;
        }
      }
    }

    for (int i = 0; i < nbins; i++) {
      if (population[igrain][i] > 0)
        autoCor[igrain][i] /= (double)population[igrain][i];
    }

    double slope = 0.0;
    for (int i = 0; i < nbins; i++) {
      double dist = i * dist_max / nbins;
      if (i > 0 && i < (nbins - 1))
        slope = 0.5 * nbins * (autoCor[igrain][i + 1] - autoCor[igrain][i - 1]) / dist_max;
      else if (i == 0)
        slope = nbins * (autoCor[igrain][1] - autoCor[igrain][0]) / dist_max;
      else
        slope = nbins * (autoCor[igrain][nbins - 1] - autoCor[igrain][nbins - 2]) / dist_max;
      file << dist << " " << autoCor[igrain][i] << " " << slope << " " << population[igrain][i] << '\n';
    }
    file << '\n';

  } // for loop igrain

  ofstream fileMeanDev("autocorMeanDev.txt");
  for (int i = 0; i < nbins; i++) {
    double Mean = 0.0;
    double fact = 1.0 / (double)num_grains;
    for (int igrain = 0; igrain < num_grains; igrain++) {
      Mean += fact * autoCor[igrain][i];
    }
    double Dev = 0.0;
    if (num_grains > 1)
      fact = 1.0 / ((double)num_grains - 1);
    for (int igrain = 0; igrain < num_grains; igrain++) {
      Dev += fact * (autoCor[igrain][i] - Mean) * (autoCor[igrain][i] - Mean);
    }
    double dist = i * dist_max / nbins;
    fileMeanDev << dist << " " << Mean << " " << sqrt(Dev) << endl;
  }
}

void pattern_fft() {
  /*
  cout << "PATTERN FOURIER TRANSFORM" << endl;
  read_image(0, iref, true); // Read image

  vector<vector<uint16_t>> patch;
  size_t lenSide = 2 * halfPatternQual + 1;
  patch.resize(lenSide);
  for (size_t l = 0; l < lenSide; l++) {
    patch[l].resize(lenSide);
  }

  for (int igrain = 0; igrain < num_grains; igrain++) {

    cout << "Patch/Mask " << igrain << endl;

    int x0 = grain[igrain].refcoord_xpix;
    int y0 = grain[igrain].refcoord_ypix;
    int xpatch = 0;
    int ypatch = 0;
    for (int i = -halfPatternQual; i <= halfPatternQual; i++) {
      ypatch = 0;
      for (int j = -halfPatternQual; j <= halfPatternQual; j++) {
        patch[xpatch][ypatch] = image[0][x0 + i][y0 + j];
        ypatch++;
      }
      xpatch++;
    }

    thumbnail patchThumb(lenSide, lenSide);
    for (size_t i = 0; i < lenSide; i++) {
      for (size_t j = 0; j < lenSide; j++) {
        patchThumb.imThumb[i][j] = patch[i][j];
      }
    }

    char fname[256];
    sprintf(fname, "patch%d.tif", igrain);
    patchThumb.writeTiff(fname, 0.0);

    // https://humbert-florent.developpez.com/algorithmique/traitement/fftw/
    size_t Npix = lenSide * lenSide;
    fftw_complex *spatial_repr = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Npix);
    fftw_complex *frequency_repr = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Npix);
    size_t offset = 0;
    for (size_t j = 0; j < lenSide; j++) {
      for (size_t i = 0; i < lenSide; i++) {
        spatial_repr[offset][0] = patch[i][j];
        spatial_repr[offset][1] = 0.0f;
        offset++;
      }
    }
    fftw_plan plan;
    plan = fftw_plan_dft_2d((int)lenSide, (int)lenSide, spatial_repr, frequency_repr, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);

    vector<float> reOut(Npix);
    vector<float> imOut(Npix);
    int x, y;
    for (int j = 0; j < (int)lenSide; j++) {
      for (int i = 0; i < (int)lenSide; i++) {
        // on recentre l'image
        x = i;
        y = j;
        if (i < halfPatternQual && j < halfPatternQual) {
          x = i + halfPatternQual;
          y = j + halfPatternQual;
        }
        if (i >= halfPatternQual && j < halfPatternQual) {
          x = i - halfPatternQual;
          y = j + halfPatternQual;
        }
        if (i < halfPatternQual && j >= halfPatternQual) {
          x = i + halfPatternQual;
          y = j - halfPatternQual;
        }
        if (i >= halfPatternQual && j >= halfPatternQual) {
          x = i - halfPatternQual;
          y = j - halfPatternQual;
        }
        reOut[y * lenSide + x] = frequency_repr[j * lenSide + i][0];
        imOut[y * lenSide + x] = frequency_repr[j * lenSide + i][1];
      }
    }

    vector<float> module(Npix);
    for (size_t i = 0; i < Npix; i++) {
      module[i] = sqrt(reOut[i] * reOut[i] + imOut[i] * imOut[i]);
    }

    float max = module[0];
    for (size_t i = 0; i < Npix / 2 - 1; i++) {
      if (max < module[i]) {
        max = module[i];
      }
    }

    if (max == 0)
      max = 1;

    for (size_t i = 0; i < Npix; i++) {
      module[i] = (module[i] / max) * UINT16_MAX;
      if (module[i] > UINT16_MAX)
        module[i] = UINT16_MAX;
    }

    thumbnail DFTThumb(lenSide, lenSide);
    for (size_t x = 0; x < lenSide; x++) {
      for (size_t y = 0; y < lenSide; y++) {
        DFTThumb.imThumb[x][y] = UINT16_MAX - module[y * lenSide + x];
      }
    }
    sprintf(fname, "FFT%d.tif", igrain);
    DFTThumb.writeTiff(fname, 0.0);

    sprintf(fname, "FFTProfile%d.txt", igrain);
    ofstream prof(fname);
    offset = halfPatternQual * lenSide;
    for (size_t x = 0; x < lenSide; x++) {
      prof << ((int)x - (int)halfPatternQual) << " " << module[offset + x] << endl;
    }

    // on détruit les objets
    fftw_destroy_plan(plan);
    fftw_free(spatial_repr);
    fftw_free(frequency_repr);

  } // loop over patches
  */
}