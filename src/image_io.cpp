
// ====================================================================
// _____ IMAGE MANAGEMENT _____
// ====================================================================

// Function that uses libtiff to read an image, and to compute and store gray levels in a c-style table.
void read_image(int i, const char *name, bool first_time) {
  double tbeg = get_time();
  fprintf(stdout, "Read image named %s... \n", name);

  // http://41j.com/blog/2011/10/simple-libtiff-example/
  TIFF *tif = TIFFOpen(name, "r");

  if (!tif) {
    cerr << "Cannot read tiff file named '" << name << "'" << endl;
    exit(0);
  }
  uint32_t w, h;
  size_t npixels;
  uint32_t *raster;
  
  
  char *infobuf;
  if (TIFFGetField(tif, TIFFTAG_DATETIME, &infobuf)) imageData[i].dateTime = std::string(infobuf);
  else imageData[i].dateTime = "dateTime unknown";
  
  imageData[i].iso_speed = 0.0;
  imageData[i].shutter = 0.0;
  imageData[i].aperture = 0.0;
  imageData[i].focal_len = 0.0;
  imageData[i].shot_order = 0; 

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  cout << w << "x" << h << endl;
  npixels = w * h;

  raster = (uint32_t *)_TIFFmalloc(npixels * sizeof(uint32_t));
  if (raster != NULL) {
    if (!TIFFReadRGBAImage(tif, w, h, raster, 0)) {
      cerr << "Cannot read tiff file named '" << name << "'" << endl;
      exit(0);
    }
  }

  if (first_time) {
    dimx = (int)w;
    dimy = (int)h;

    // Reserve memory for image
    if (!image.empty())
      image.clear();
    image.resize(2);
    for (size_t i = 0; i < 2; i++) {
      image[i].resize(dimx);
      for (int ix = 0; ix < dimx; ix++) {
        image[i][ix].resize(dimy);
      }
    }
  }

  double r, g, b;
  double fact = 1.0 / 255.0;
  for (int row = 0; row < dimy; ++row) {
    int shift = row * w;
    for (int col = 0; col < dimx; ++col) {
      int p = shift + col;
      r = TIFFGetR(raster[p]) * fact;
      g = TIFFGetG(raster[p]) * fact;
      b = TIFFGetB(raster[p]) * fact;
      if (r < 0.0) { r = 0.0; std::cerr << "@read_image, r < 0.0\n"; }
      if (r > 1.0) { r = 1.0; std::cerr << "@read_image, r > 1.0\n"; }
      if (g < 0.0) { g = 0.0; std::cerr << "@read_image, g < 0.0\n"; }
      if (g > 1.0) { g = 1.0; std::cerr << "@read_image, g > 1.0\n"; }
      if (b < 0.0) { b = 0.0; std::cerr << "@read_image, b < 0.0\n"; }
      if (b > 1.0) { b = 1.0; std::cerr << "@read_image, b > 1.0\n"; }
      // r, g and b are now values in the range [0, 1]

      // Store in a c-style table, and convert in graylevel
      // Remember that in tiff format the coordinate (0, 0) is in lower-left corner.
      image[i][col][h - row - 1] = (uint16_t)(65535.0 * (0.299 * r + 0.587 * g + 0.114 * b));
    }
  }

  _TIFFfree(raster);
  TIFFClose(tif);

  if (first_time) {
    int j = 1 - i;
    for (int y = 0; y < dimy; ++y) {
      for (int x = 0; x < dimx; ++x) {
        image[j][x][y] = image[i][x][y];
      }
    }
  }

  fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);
}

// il faut voir ici http://www.libraw.org/node/555 pour voir comment proceder
// ici aussi :
// http://stackoverflow.com/questions/22355491/libraw-is-making-my-images-too-bright-compared-to-nikons-own-converter
void read_raw_image(int i, const char *name, bool first_time) {
#ifndef CYGWIN
  LibRaw iProcessor;
  int IO_error = iProcessor.open_file(name);
  if (IO_error != 0) {
    cout << endl;
    cout << "Sorry but your image named " << name << " does not exist..." << endl;
    cout << "   ..PROGRAM STOPPED.." << endl;
    exit(1);
  }

  double tbeg = get_time();
  fprintf(stdout, "Read and demosaicing image named %s (DemosaicModel = %d)... \n", name, DemosaicModel);

  if (first_time) {
    dimx = iProcessor.imgdata.sizes.width;
    dimy = iProcessor.imgdata.sizes.height;
    cout << "Size = " << dimx << "x" << dimy << endl;

    // Reserve memory for image
    if (!image.empty())
      image.clear();
    image.resize(2);
    for (size_t im = 0; im < 2; im++) {
      image[im].resize(dimx);
      for (int ix = 0; ix < dimx; ix++) {
        image[im][ix].resize(dimy);
      }
    }
  }
  
  // Get data about the image shot
  imageData[i].iso_speed = iProcessor.imgdata.other.iso_speed;
  imageData[i].shutter = iProcessor.imgdata.other.shutter;
  imageData[i].aperture = iProcessor.imgdata.other.aperture;
  imageData[i].focal_len = iProcessor.imgdata.other.focal_len;
  imageData[i].shot_order = iProcessor.imgdata.other.shot_order;
  imageData[i].dateTime = timestamp2string(iProcessor.imgdata.other.timestamp);
  
  iProcessor.unpack();
  uint16_t MinGray = 65535, MaxGray = 0;
  double fact = 1.0 / (double)(iProcessor.imgdata.color.maximum);
  //std::cout << "iProcessor.imgdata.color.maximum = " << iProcessor.imgdata.color.maximum << '\n';
  //std::cout << "fact = " << fact << '\n';

  if (DemosaicModel >= 0) { // That are those of libRaw
    iProcessor.imgdata.params.user_qual = DemosaicModel;
    iProcessor.dcraw_process();
    double r, g, b;
    int yshift = 0;
    for (int y = 0; y < dimy; ++y) {
      yshift = y * dimx;
      for (int x = 0; x < dimx; ++x) {
        // According to http://www.libraw.org/node/1995, we can ignore g2 level
        int pos = yshift + x;
        r = iProcessor.imgdata.image[pos][0] * fact;
        g = iProcessor.imgdata.image[pos][1] * fact;
        b = iProcessor.imgdata.image[pos][2] * fact;

        // Store in a c-style table, and convert in graylevel
        image[i][x][y] = (uint16_t)(iProcessor.imgdata.color.maximum * (0.299 * r + 0.587 * g + 0.114 * b));

        if (image[i][x][y] > MaxGray)
          MaxGray = image[i][x][y];
        if (image[i][x][y] < MinGray)
          MinGray = image[i][x][y];
      }
    }
  } else if (DemosaicModel == -1) {

    int width = iProcessor.imgdata.sizes.raw_width;
    int yoffset = iProcessor.imgdata.sizes.top_margin;
    int xoffset = iProcessor.imgdata.sizes.left_margin;

    double deraster;
    for (int x = 0; x < dimx; ++x) {
      for (int y = 0; y < dimy; ++y) {
        // 11
        // 11
        deraster = 0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y)]);
        deraster += 0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y)]);
        deraster +=
            0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y + 1)]);
        deraster += 0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y + 1)]);

        image[i][x][y] = (unsigned int)floor(deraster);
        if (image[i][x][y] > MaxGray)
          MaxGray = image[i][x][y];
        if (image[i][x][y] < MinGray)
          MinGray = image[i][x][y];
      }
    }
  } else if (DemosaicModel == -2) {

    int width = iProcessor.imgdata.sizes.raw_width;
    int yoffset = iProcessor.imgdata.sizes.top_margin;
    int xoffset = iProcessor.imgdata.sizes.left_margin;

    double deraster;
    for (int x = 0; x < dimx; ++x) {
      for (int y = 0; y < dimy; ++y) {
        // 121
        // 242
        // 121
        deraster =
            0.0625 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x - 1) + width * (yoffset + y - 1)]);
        deraster += 0.125 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y - 1)]);
        deraster +=
            0.0625 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y - 1)]);

        deraster += 0.125 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x - 1) + width * (yoffset + y)]);
        deraster += 0.25 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y)]);
        deraster += 0.125 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y)]);

        deraster +=
            0.0625 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x - 1) + width * (yoffset + y + 1)]);
        deraster += 0.125 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y + 1)]);
        deraster +=
            0.0625 * (double)(iProcessor.imgdata.rawdata.raw_image[(xoffset + x + 1) + width * (yoffset + y + 1)]);

        image[i][x][y] = (unsigned int)floor(deraster);
        if (image[i][x][y] > MaxGray)
          MaxGray = image[i][x][y];
        if (image[i][x][y] < MinGray)
          MinGray = image[i][x][y];
      }
    }
  } else if (DemosaicModel == -20) { // for achromatic raw images (no demosaicing)

    int width = iProcessor.imgdata.sizes.raw_width;
    int yoffset = iProcessor.imgdata.sizes.top_margin;
    int xoffset = iProcessor.imgdata.sizes.left_margin;

    for (int x = 0; x < dimx; ++x) {
      for (int y = 0; y < dimy; ++y) {
        image[i][x][y] = iProcessor.imgdata.rawdata.raw_image[(xoffset + x) + width * (yoffset + y)];
        if (image[i][x][y] > MaxGray)
          MaxGray = image[i][x][y];
        if (image[i][x][y] < MinGray)
          MinGray = image[i][x][y];
      }
    }
    
  }

  iProcessor.recycle();
  
  std::cout << "MinGray = " << MinGray << '\n';
  std::cout << "MaxGray = " << MaxGray << '\n';
  
  // rescale according to MaxGray and MinGray. 
  // It should be ok since we use ZNCC (which is not sensitive to a scale factor)
  // The default behaviour is to NOT rescale the gray levels
  if (rescaleGrayLevels) {
    std::cout << "rescaling !!!\n";
    double fact = 65535.0 / (double)(MaxGray - MinGray);
    for (int y = 0; y < dimy; ++y) {
      for (int x = 0; x < dimx; ++x) {
        image[i][x][y] = (uint16_t)floor((image[i][x][y] - MinGray) * fact);
      }
    }
  }

  if (first_time) {
    int j = 1 - i;
    for (int y = 0; y < dimy; ++y) {
      for (int x = 0; x < dimx; ++x) {
        image[j][x][y] = image[i][x][y];
      }
    }
  }

  fprintf(stdout, "[DONE in %f seconds, user_qual = %d]\n", get_time() - tbeg, iProcessor.imgdata.params.user_qual);

#if 0
	int x0_d = 10;
	int y0_d = 10;
	int W_d = 2000;
	int H_d = 2000;
	double rr = 255./65535.;
	ofstream file("demosaicing.pgm", ios::binary);
	file << "P5\n";
	file << W_d << " " << H_d << endl;
	file << "255" << endl;

	for (int y = 0 ; y < H_d ; y++) {
		for (int x = 0 ; x < W_d ; x++) {
			file << (int)(rr*image[i][x + x0_d][y + y0_d]) << endl;
		}
	}

	exit(0);
#endif
#endif
}

// This function has to be rewritten for a pgm (or tiff) output
void undistor_image(const char * /*name_from*/, const char * /*name_to*/) {  
  read_image(im_index_ref, grid_image_name.c_str(), true);
  
  double xu, yu;
	for (int y = 0 ; y < dimy ; ++y) {
		for (int x = 0 ; x < dimx ; ++x) {
			undistor(&disto_parameters[0], x, y, xu, yu);
			xu = nearest(xu);
			yu = nearest(yu);
			if (xu < 0. || xu >= dimx || yu < 0. || yu >= dimy) continue;
      image[im_index_current][(int)xu][(int)yu] = image[im_index_ref][x][y];
		}
	}
  
  create_image_Netbpm(0);
}

// This is the main function to read an image.
// The image can be RAW or TIFF so that the library libraw or libtiff is used.
// The image format is not recognized, but the user needs to set RawImage=1 to say that RAW format has to be used
void read_image(int i, int num, bool first_time) {
  char name[256];
  sprintf(name, image_name, num);
  if (RawImages)
    read_raw_image(i, name, first_time);
  else
    read_image(i, name, first_time);
}

void draw_circle(vector<vector<uint16_t>> &imThumb, int x_centre, int y_centre, int radius, uint16_t col) {
  int x = 0, y = radius, m = 5 - 4 * radius;
  while (x <= y) { /// @see http://fr.wikipedia.org/wiki/Algorithme_de_tracé_d'arc_de_cercle_de_Bresenham
    imThumb[x + x_centre][y + y_centre] = col;
    imThumb[y + x_centre][x + y_centre] = col;
    imThumb[-x + x_centre][y + y_centre] = col;
    imThumb[-y + x_centre][x + y_centre] = col;
    imThumb[x + x_centre][-y + y_centre] = col;
    imThumb[y + x_centre][-x + y_centre] = col;
    imThumb[-x + x_centre][-y + y_centre] = col;
    imThumb[-y + x_centre][-x + y_centre] = col;
    if (m > 0) {
      y = y - 1;
      m = m - 8 * y;
    }
    x = x + 1;
    m = m + 8 * x + 4;
  }
}

// An alternative solution to make helper snapshots if Magick++ is not used
// image_div is a global parameter that is use to divide the side sizes of the image
void create_image_Netbpm(int num) {
  // Create an resized image
  int imW = dimx / image_div; // integer division
  int imH = dimy / image_div;
  cout << "Thumb size: " << imW << "x" << imH << endl;

  vector<vector<uint16_t>> imThumb;
  imThumb.resize(imW);
  for (int i = 0; i < imW; i++) {
    imThumb[i].resize(imH);
  }

  // Resize
  int x0 = 0, y0 = 0;
  uint16_t mean_grey;
  double fact;
  int ix = 0, iy = 0;

  uint16_t max_grey = 0;
  uint16_t min_grey = 65535; // suppose 16-bits
  for (int xm = image_div; xm < dimx; xm += image_div) {
    y0 = 0;
    iy = 0;
    for (int ym = image_div; ym < dimy; ym += image_div) {
      mean_grey = 0;
      fact = 1.0 / ((xm - x0) * (ym - y0));
      for (int x = x0; x < xm; x++) {
        for (int y = y0; y < ym; y++) {
          mean_grey += nearest(fact * image[im_index_current][x][y]);
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
  cout << "min_grey = " << min_grey << endl;
  cout << "max_grey = " << max_grey << endl;
  double rescaleFactor = 1.0;
  if ((max_grey - min_grey) > 0)
    rescaleFactor = 65535.0 / (max_grey - min_grey);
  for (int y = 0; y < imH; y++) {
    for (int x = 0; x < imW; x++) {
      imThumb[x][y] = (uint16_t)nearest((imThumb[x][y] - min_grey) * rescaleFactor);
    }
  }

  // Draw disks
  uint16_t disk_grey = 65535;
  for (int i = 0; i < num_grains; ++i) {
    int x_centre = (int)(grain[i].refcoord_xpix + grain[i].dx);
    int y_centre = (int)(grain[i].refcoord_ypix + grain[i].dy);
    int rayon = (int)grain[i].radius_pix;

    x_centre /= image_div;
    y_centre /= image_div;
    rayon /= image_div;

    draw_circle(imThumb, x_centre, y_centre, rayon, disk_grey);

    if (grain[i].NCC_rescue <= NCC_min) { // Un bricolage pour identifier le grain fautife
      int x = 0, y = rayon, m = 5 - 4 * rayon;
      while (x <= y) {
        for (int i = 0; i < x; i++) {
          imThumb[i + x_centre][y + y_centre] = disk_grey;
          imThumb[y + x_centre][i + y_centre] = disk_grey;
          imThumb[-i + x_centre][y + y_centre] = disk_grey;
          imThumb[-y + x_centre][i + y_centre] = disk_grey;
          imThumb[i + x_centre][-y + y_centre] = disk_grey;
          imThumb[y + x_centre][-i + y_centre] = disk_grey;
          imThumb[-i + x_centre][-y + y_centre] = disk_grey;
          imThumb[-y + x_centre][-i + y_centre] = disk_grey;
        }

        if (m > 0) {
          y = y - 1;
          m = m - 8 * y;
        }
        x = x + 1;
        m = m + 8 * x + 4;
      }
    }
  }

  // Write the file
  cout << "Create image " << num << " for check ... " << flush;
  char name[256];
  sprintf(name, "check_%d.pgm", num);
  ofstream file(name, ios::binary);
  file << "P5\n";
  file << imW << " " << imH << endl;
  file << "255" << endl;

  double conv = 255.0 / 65535.0;

  for (int y = 0; y < imH; y++) {
    for (int x = 0; x < imW; x++) {
      unsigned char v = nearest(conv * (imThumb[x][y]));
      file.write((const char *)&v, sizeof(unsigned char));
    }
  }

  cout << "done." << endl;
}

void create_image(int num) {
  create_image_Netbpm(num);
}

