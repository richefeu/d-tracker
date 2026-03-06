// -----------------------------------
// TRACKER 2.x, December 2011
// Gael.Combe@ujf-grenoble.fr
// Vincent.Richefeu@ujf-grenoble.fr
// Lab 3SR, Grenoble
// -----------------------------------

#define TRACKER_VERSION "0.6.0"

#include "tracker.hpp"

void header() {
  // http://www.patorjk.com/software/taag
  // Lean
  cout << endl;
  cout << "_/_/_/_/_/  _/_/_/      _/_/      _/_/_/  _/    _/  _/_/_/_/  _/_/_/ " << endl;
  cout << "   _/      _/    _/  _/    _/  _/        _/  _/    _/        _/    _/" << endl;
  cout << "  _/      _/_/_/    _/_/_/_/  _/        _/_/      _/_/_/    _/_/_/   " << endl;
  cout << " _/      _/    _/  _/    _/  _/        _/  _/    _/        _/    _/  " << endl;
  cout << "_/      _/    _/  _/    _/    _/_/_/  _/    _/  _/_/_/_/  _/    _/   " << endl;
  cout << endl;
  cout << " by Vincent.Richefeu@3sr-grenoble.fr AND Gael.Combe@3sr-grenoble.fr  " << endl;
  cout << endl << endl;
  cout << "_____ Version number: " << TRACKER_VERSION << endl;

#if defined(_OPENMP)
  cout << "_____ Accelerator: OPENMP" << endl;
#else
  cout << "_____ Accelerator: None" << endl;
#endif

#if defined(SVN_REV)
  cout << "_____ SVN revision number: " << SVN_REV << endl;
#endif
  cout << endl << endl;
}

double get_time() {
#if defined(_OPENMP)
  return omp_get_wtime();
#else
  return (double)std::clock() / (double)CLOCKS_PER_SEC;
#endif
}

void dialog() {
  header();
  cout << endl;
  cout << " [ Enter Interactive Mode ]" << endl;
  cout << endl;
  cout << "Type 'help' for help" << endl;

  bool loaded = false;

  string token;
  cout << "tracker (interactive mode) > ";
  cin >> token;
  for (;;) {
    if (token == "help") {
      cout << "Type one of the following keywords:" << endl;
      cout << "  read_image  ..................... Read an image by using libtiff" << endl;
      cout << "  read_raw_image  ................. Read an image by using libraw" << endl;
      cout << "  DemosaicModel  .................. Set the value of DemosaicModel" << endl;
      cout << "  minmax  ......................... Show min and max values in the loaded image" << endl;
      cout << "  histo ........................... Compute the histogram of the loaded image" << endl;
      cout << "  image_check ..................... Create an image as it is 'seen' by tracker" << endl;
      // cout << "  draw_grains ..................... Superimpose grains over the thumb" << endl;
    } else if (token == "test") {
      im_index_current = 0;
      // read_raw_image(im_index_current, "../examples/ImageData/Disto_1.IIQ", true);
      read_image(im_index_current, "../examples/ImageData/TEST02.tiff", true);
      read_grains("../examples/ParticleImageTracking/dic_out_2.txt", true);
      image_div = 1;
      thumbnail tt(image, im_index_current, image_div);
      ColorTable ct;
      ct.setMinMax(0.0, 1.0);
      ct.SavePpm("ColorScale.ppm");
      colorRGBA cRGBA;
      for (int i = 0; i < num_grains; ++i) {
        int x_centre = (int)(grain[i].refcoord_xpix + grain[i].dx);
        int y_centre = (int)(grain[i].refcoord_ypix + grain[i].dy);
        int rayon = (int)grain[i].radius_pix;

        x_centre /= image_div;
        y_centre /= image_div;
        rayon /= image_div;

        ct.getRGB(grain[i].NCC, &cRGBA);
        RGB col;
        col.r = cRGBA.r;
        col.g = cRGBA.g;
        col.b = cRGBA.b;
        tt.draw_disk(x_centre, y_centre, rayon, col);
      }
      tt.writePpm("visu.ppm", 0.3);
      tt.writeTiff("visu.tiff", 0.3);
    } else if (token == "test_raw") {
      im_index_current = 0;
      read_raw_image(im_index_current, "../examples/ImageData/CaptureOneSession2687.IIQ", true);
      image_div = 4;
      thumbnail tt(image, im_index_current, image_div);
      // tt.writePpm("visu.ppm", 0.3);
      tt.writeTiff("visuRaw.tiff", 0.0);
    } else if (token == "read_dic") {
      cout << "file path: ";
      string path;
      cin >> path;
      if (!grain.empty())
        grain.clear();
      read_grains(path.c_str(), true);
    } else if (token == "read_image") {
      cout << "file path: ";
      string path;
      cin >> path;
      im_index_current = 0;
      read_image(im_index_current, path.c_str(), true);
      if (dimx != 0)
        loaded = true;
    } else if (token == "read_raw_image") {
      cout << "file path: ";
      string path;
      cin >> path;
      im_index_current = 0;
      read_raw_image(im_index_current, path.c_str(), true);
      if (dimx != 0)
        loaded = true;
    } else if (token == "DemosaicModel") {
      cout << "value: ";
      cin >> DemosaicModel;
    } else if (token == "image_check") {
      if (!loaded) {
        cout << "Read an image first" << endl;
        token = "";
        continue;
      }
      im_index_current = 0;
      image_div = 1;
      thumbnail tt(image, im_index_current, image_div);
      tt.writeTiff("image_check.tiff", 0.0);
    } else if (token == "minmax") {
      if (!loaded) {
        cout << "Read an image first" << endl;
        token = "";
        continue;
      }
      uint16_t minVal = 65535; // suppose 16-bit depth
      uint16_t maxVal = 0;
      for (int ix = 0; ix < dimx; ix++) {
        for (int iy = 0; iy < dimy; iy++) {
          if (image[0][ix][iy] > maxVal)
            maxVal = image[0][ix][iy];
          if (image[0][ix][iy] < minVal)
            minVal = image[0][ix][iy];
        }
      }
      cout << "minVal = " << minVal << endl;
      cout << "maxVal = " << maxVal << endl;
    } else if (token == "histo") {
      if (!loaded) {
        cout << "Read an image first" << endl;
        token = "";
        continue;
      }

      vector<size_t> hist(65535);
      for (size_t i = 0; i < 65535; i++) {
        hist[i] = 0;
      }

      for (int ix = 0; ix < dimx; ix++) {
        for (int iy = 0; iy < dimy; iy++) {
          hist[image[0][ix][iy]]++;
        }
      }

      ofstream file("histo.txt");
      for (size_t i = 0; i < 65535; i++) {
        file << hist[i] << endl;
      }
      cout << "Histogram saved in 'hist.txt'" << endl;
    } else if (token == "q" || token == "quit")
      break;

    cout << "tracker> ";
    cin >> token;
  }
  cout << " [ Exit Interactive Mode ]" << endl;
  exit(EXIT_SUCCESS);
}

/************************************************************************************************/
/*  				 PARTICLE TRACKING                                                                  */
/************************************************************************************************/

// The particle tracking procedure (default procedure)
void particle_tracking() {
  read_image(im_index_ref, iref, true); // Read reference image

  if (make_images)
    create_image(iref);

  double tbeg;
  int num_image;
  int igrain;
  int irescue;

  if (use_neighbour_list) {
    LOG("Build neighbour list ... " << flush);
    tbeg = get_time();
    for (int igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].neighbour.clear();
      find_neighbours(igrain);
    }
    LOG("[DONE in " << get_time() - tbeg << " seconds]" << endl);
  }

  for (num_image = ibeg; num_image <= iend; num_image += iinc) {

    fprintf(stdout, "\n____ Correlations from image %d towards image %d (%d/%d)\n", iref, num_image,
            num_image - ibeg + 1, iend - ibeg + 1);
    read_image(im_index_current, num_image);

    // Re-build the neighbour list periodically
    if (use_neighbour_list && num_image % period_rebuild_neighbour_list == 0) {
      fprintf(stdout, "Rebuilding the neighbour list ... ");
      fflush(stdout);
      tbeg = get_time();
      for (int igrain = 0; igrain < num_grains; igrain++) {
        find_neighbours(igrain);
      }

      fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);
    }

    // Precomputations
    if (require_precomputations) {
      fprintf(stdout, "Precomputations ... ");
      fflush(stdout);
      tbeg = get_time();
      do_precomputations();
      fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);
    }

    // non sub-pixel correlation (with rotations).
    // The very first attempts
    num_to_be_rescued = num_to_be_super_rescued = 0;
    fprintf(stdout, "Follow %d grains ... \n", num_grains);
    fflush(stdout);

    if (rescue_level >= 0) {
      tbeg = get_time();

      progress = 0;
#pragma omp parallel for schedule(dynamic)
      for (igrain = 0; igrain < num_grains; igrain++) {
        grain[igrain].reset();  // reset les NCC
        grain[igrain].backup(); // sauvegarde les dx, dy et drot
        if (grain[igrain].masked)
          continue;
        follow_pattern_pixel(igrain);
#pragma omp critical
        { loadbar(++progress, grain.size()); }
      }

      cerr << endl;

      fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);
    }

    // According to a minimum value of NCC, a rescue is attempted
    if (rescue_level >= 1) {
      if (num_to_be_rescued > 0)
        fprintf(stdout, "Try to rescue %d grains\n", num_to_be_rescued);

      tbeg = get_time();

      progress = 0;
#pragma omp parallel for private(igrain) schedule(dynamic)
      for (irescue = 0; irescue < num_to_be_rescued; irescue++) {
        igrain = to_be_rescued[irescue];
        follow_pattern_rescue_pixel(igrain);
#pragma omp critical
        { loadbar(++progress, num_to_be_rescued); }
      }
      cerr << endl;

      // Display
      for (irescue = 0; irescue < num_to_be_rescued; irescue++) {
        igrain = to_be_rescued[irescue];
        fprintf(stdout, "\nGrain %6d \t NCC = %6.4f  ->  ", igrain, grain[igrain].NCC);
        fprintf(stdout, "%6.4f \t diff = %6.4f \t ", grain[igrain].NCC_rescue,
                grain[igrain].NCC_rescue - grain[igrain].NCC);
        if (grain[igrain].NCC_rescue > NCC_min)
          fprintf(stdout, "[\033[32mOK\033[0m]\n");
        else
          fprintf(stdout, "[\033[31mFAIL\033[0m]\n");
      }

      if (num_to_be_rescued > 0)
        fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);
    }

    // According to a minimum value of NCC, a last super_rescue is attempted
    if (rescue_level >= 2) {
      if (num_to_be_super_rescued > 0)
        fprintf(stdout, "Try to super_rescue %d grains\n", num_to_be_super_rescued);

      tbeg = get_time();

      progress = 0;
#pragma omp parallel for private(igrain) schedule(dynamic)
      for (irescue = 0; irescue < num_to_be_super_rescued; irescue++) {
        igrain = to_be_super_rescued[irescue];
        follow_pattern_super_rescue_pixel(igrain);
#pragma omp critical
        { loadbar(++progress, num_to_be_rescued); }
      }
      cerr << endl;

      // Display
      for (irescue = 0; irescue < num_to_be_super_rescued; irescue++) {
        igrain = to_be_super_rescued[irescue];
        fprintf(stdout, "Grain %6d \t NCC = %6.4f  ->  ", igrain, grain[igrain].NCC);
        fprintf(stdout, "%6.4f \t diff = %6.4f \t ", grain[igrain].NCC_rescue,
                grain[igrain].NCC_rescue - grain[igrain].NCC);
        if (grain[igrain].NCC_rescue > NCC_min_super)
          fprintf(stdout, "[\033[32mOK\033[0m]\n");
        else {
          fprintf(stdout, "[\033[31mFAIL\033[0m]\n");
          fprintf(stdout, " Please, fix the pb manually\n");
        }
      }

      if (num_to_be_super_rescued > 0)
        fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);
    }

    // Sub-pixel
    if (subpixel) {
      fprintf(stdout, "Sub-pixel resolution ... \n");
      fflush(stdout);

      tbeg = get_time();

      progress = 0;
#pragma omp parallel for schedule(dynamic)
      for (igrain = 0; igrain < num_grains; igrain++) {
        if (grain[igrain].masked)
          continue;
        follow_pattern_subpixel_xyR(igrain);
#pragma omp critical
        { loadbar(++progress, num_grains); }
      }
      cerr << endl;
      fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);
    }

    for (igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].upix = grain[igrain].dx - grain[igrain].dx_prev;
      grain[igrain].vpix = grain[igrain].dy - grain[igrain].dy_prev;
      grain[igrain].rot_inc = grain[igrain].drot - grain[igrain].drot_prev;
    }

    // Save
    save_grains(num_image);

    if (make_images)
      create_image(num_image);

    // Change reference (Not yet tested!!!)
    if (num_image - iref >= iraz) { // equal in fact!
      cout << endl;
      cout << "- - - - - - - - - - - - - - - - - - - - - -" << endl << endl;
      cout << "Changing reference image from " << iref << " to " << num_image << endl << endl;
      cout << "- - - - - - - - - - - - - - - - - - - - - -" << endl;

      int tmp = im_index_current;
      im_index_current = im_index_ref;
      im_index_ref = tmp; // Swap the first index for the array 'image'
      double actual_pos;
      for (igrain = 0; igrain < num_grains; igrain++) {
        actual_pos = (double)(grain[igrain].refcoord_xpix) + grain[igrain].dx;
        grain[igrain].refcoord_xpix = nearest(actual_pos);
        grain[igrain].dx = actual_pos - grain[igrain].refcoord_xpix;

        actual_pos = (double)(grain[igrain].refcoord_ypix) + grain[igrain].dy;
        grain[igrain].refcoord_ypix = nearest(actual_pos);
        grain[igrain].dy = actual_pos - grain[igrain].refcoord_ypix;

        grain[igrain].refrot += grain[igrain].drot;
        grain[igrain].drot = 0.0;
      }
      iref = num_image;
      require_precomputations = true; // so that precomputations are updated
    }
  }

  if (ShowSolidMotionError) {
    double angle = 0.0, xc = 0.0, yc = 0.0, xtrans = grain[0].dx, ytrans = grain[0].dy;
    double mean = 0.0, stddev = 0.0;
    SolidMotionError(angle, xc, yc, xtrans, ytrans, mean, stddev);
    cout << "@SolidMotionError, meanError = " << mean << ", stddevError = " << stddev << endl;
    cout << "angle  = " << angle << endl;
    cout << "xc     = " << xc << endl;
    cout << "yc     = " << yc << endl;
    cout << "xtrans = " << xtrans << endl;
    cout << "ytrans = " << ytrans << endl;
  }
}

// Une version de particle tracking qui utilise une seconde strategie (en test pour le moment)
// L'idée est de ne faire que du subpix !!!!
//
void particle_tracking_v2_subpix() {
  // TODO !!!!!!!!!
}

void rotate_pixel_pattern(int igrain, int i, double c, double s, int *xpixel, int *ypixel) {
  // Attention y est oriente vers le bas, c'est pourquoi le sens trigo est inverse
  double xpix = c * grain[igrain].pattern[i].dx + s * grain[igrain].pattern[i].dy;
  double ypix = -s * grain[igrain].pattern[i].dx + c * grain[igrain].pattern[i].dy;
  *xpixel = nearest(xpix);
  *ypixel = nearest(ypix);
}

// Pre-computation of mean0, C0C0 and interpolated_zone for each grain
void do_precomputations() {
  double c0, s0;
  int xpixel0, ypixel0;
  double diffC0;

  for (int igrain = 0; igrain < num_grains; igrain++) {
    int refcoordx = grain[igrain].refcoord_xpix;
    int refcoordy = grain[igrain].refcoord_ypix;
    double refrot = grain[igrain].refrot;

    grain[igrain].mean0 = 0.0;
    grain[igrain].C0C0 = 0.0;
    c0 = cos(refrot);
    s0 = sin(refrot);
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      rotate_pixel_pattern(igrain, i, c0, s0, &xpixel0, &ypixel0);
      grain[igrain].pattern0_rotated[i].dx = xpixel0;
      grain[igrain].pattern0_rotated[i].dy = ypixel0;
      xpixel0 += refcoordx;
      ypixel0 += refcoordy;
      grain[igrain].mean0 += (double)image[im_index_ref][xpixel0][ypixel0];
    }
    grain[igrain].mean0 /= (double)grain[igrain].pattern.size();
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      xpixel0 = grain[igrain].pattern0_rotated[i].dx + refcoordx;
      ypixel0 = grain[igrain].pattern0_rotated[i].dy + refcoordy;
      diffC0 = (double)image[im_index_ref][xpixel0][ypixel0] - grain[igrain].mean0;
      grain[igrain].C0C0 += diffC0 * diffC0;
    }
  } // Loop over igrain

#if 0
	// Utile pour faire un graphe avec des valeurs interpolees
	ofstream bidon("bidon1.txt");
	for (double x = grain[0].refcoord_xpix - 3 ; x <= grain[0].refcoord_xpix + 3 ; x += 0.1) {
		for (double y = grain[0].refcoord_ypix - 3 ; y <= grain[0].refcoord_ypix + 3 ; y += 0.1) {
			bidon << x << " " << y << " " << IMAGE_INTERPOLATOR->getValue(im_index_current, x, y) << endl;
		}
		bidon << endl;
	}

	ofstream bidon2("bidon2.txt");

	for (int x = grain[0].refcoord_xpix - 3 ; x <= grain[0].refcoord_xpix + 3 ; x++) {
		for (int y = grain[0].refcoord_ypix - 3 ; y <= grain[0].refcoord_ypix + 3 ; y++) {
			bidon2 << x << " " << y << " " << image[im_index_ref][x][y] << endl;
		}
	}
	exit(0);
#endif

  require_precomputations = false;
}

// Based on Highest NCC
void follow_pattern_pixel(int igrain) {
  int refcoordx = grain[igrain].refcoord_xpix;
  int refcoordy = grain[igrain].refcoord_ypix;
  double refrot = grain[igrain].refrot;
  double c1, s1;
  int xpixel0, ypixel0;

  int grain_dx = nearest(grain[igrain].dx);
  int grain_dy = nearest(grain[igrain].dy);
  double rest_dx = grain[igrain].dx - grain_dx;
  double rest_dy = grain[igrain].dy - grain_dy;
  grain[igrain].upix = grain[igrain].vpix = grain[igrain].rot_inc = 0.0;

  double drot;    // Increment test pour l'angle de rotation
  int upix, vpix; // Increment test pour la translation

  double mean1, C0C1, C1C1;

  int xpixel1, ypixel1;
  double total_rot;

  double NCC_test, best_NCC = 0.0;
  double best_drot = 0.0;
  int best_upix = 0.0, best_vpix = 0.0;
  double diffC1;
  double rotmax = fabs(search_zone.inc_rot * search_zone.num_rot);

  for (drot = -rotmax; drot <= rotmax; drot += search_zone.inc_rot) {
    total_rot = refrot + grain[igrain].drot + drot;
    c1 = cos(total_rot);
    s1 = sin(total_rot);
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      rotate_pixel_pattern(igrain, i, c1, s1, &xpixel1, &ypixel1);
      grain[igrain].pattern1_rotated[i].dx = xpixel1;
      grain[igrain].pattern1_rotated[i].dy = ypixel1;
    }
    for (upix = -search_zone.left; upix <= search_zone.right; upix++) {
      for (vpix = -search_zone.up; vpix <= search_zone.down; vpix++) {

        mean1 = 0.0;
        C0C1 = C1C1 = 0.0;
        for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
          xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
          ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
          mean1 += (double)image[im_index_current][xpixel1][ypixel1];
        }
        mean1 /= (double)grain[igrain].pattern.size();
        for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
          xpixel0 = grain[igrain].pattern0_rotated[i].dx + refcoordx;
          ypixel0 = grain[igrain].pattern0_rotated[i].dy + refcoordy;
          xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
          ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
          diffC1 = (double)image[im_index_current][xpixel1][ypixel1] - mean1;
          C0C1 += ((double)image[im_index_ref][xpixel0][ypixel0] - grain[igrain].mean0) * diffC1;
          C1C1 += diffC1 * diffC1;
        }

        NCC_test = C0C1 / sqrt(grain[igrain].C0C0 * C1C1);

        if (best_NCC < NCC_test) {
          best_NCC = NCC_test;
          best_upix = (double)upix;
          best_vpix = (double)vpix;
          best_drot = drot;
        }

      } // Loop vpix
    }   // Loop upix
  }     // Loop rot

  grain[igrain].NCC = best_NCC;
  grain[igrain].NCC_rescue = best_NCC;

  if (best_NCC < NCC_min) {
#pragma omp critical
    { to_be_rescued[num_to_be_rescued++] = igrain; }
  } else {
    grain[igrain].upix = best_upix - rest_dx;
    grain[igrain].vpix = best_vpix - rest_dy;
    grain[igrain].rot_inc = best_drot;
    grain[igrain].dx = grain_dx + best_upix;
    grain[igrain].dy = grain_dy + best_vpix;
    grain[igrain].drot += best_drot;
  }
}

void follow_pattern_rescue_pixel(int igrain) {
  int refcoordx = grain[igrain].refcoord_xpix;
  int refcoordy = grain[igrain].refcoord_ypix;
  double refrot = grain[igrain].refrot;
  double c, s;
  int xpixel0, ypixel0;

  int grain_dx = nearest(grain[igrain].dx); // modif by gael 15/01/2012
  int grain_dy = nearest(grain[igrain].dy); // modif by gael 15/01/2012
  double rest_dx = grain[igrain].dx - grain_dx;
  double rest_dy = grain[igrain].dy - grain_dy;

  double drot;    // Increment test pour l'angle de rotation
  int upix, vpix; // Increment test pour la translation

  double mean1, C0C1, C1C1;

  int i_allowed;
  bool is_allowed;

  int xpixel1, ypixel1;
  double total_rot;

  double NCC_test, best_NCC = 0.0;
  double best_drot = 0.0;
  int best_upix = 0.0, best_vpix = 0.0;
  double rotmax = fabs(search_zone_rescue.inc_rot * search_zone_rescue.num_rot);

  // We define a list of relative coordinates used to search igrain in the image 1.
  // To do that, the list of neighbour of igrain is used
  const int nbmax = (search_zone_rescue.right + search_zone_rescue.left + 1) *
                    (search_zone_rescue.up + search_zone_rescue.down + 1);
  vector<relative_coord_type> rescue_allowed(nbmax);

  int radi_i = (int)grain[igrain].radius_pix;

  int nb_allowed = 0;

  if (use_neighbour_list) {
    for (int ix = -search_zone_rescue.left; ix <= search_zone_rescue.right; ix++) {
      for (int iy = -search_zone_rescue.up; iy <= search_zone_rescue.down; iy++) {
        is_allowed = true;
        int xcentre_i = refcoordx + grain_dx + ix;
        int ycentre_i = refcoordy + grain_dy + iy;
        for (int i_neigh = 0; i_neigh < grain[igrain].num_neighbour; i_neigh++) {
          int jgrain = grain[igrain].neighbour[i_neigh];
          double jgrain_NCC = grain[jgrain].NCC_rescue;
          if (jgrain_NCC > NCC_min) {
            int xcentre_j = grain[jgrain].refcoord_xpix + nearest(grain[jgrain].dx);
            int ycentre_j = grain[jgrain].refcoord_ypix + nearest(grain[jgrain].dy);
            int radi_j = (int)grain[jgrain].radius_pix;
            double dstx = xcentre_j - xcentre_i;
            double dsty = ycentre_j - ycentre_i;
            double dist2 = dstx * dstx + dsty * dsty;
            double sum_radii = (radi_i + radi_j) * 0.9; // only 90% of the radius declared
            double sum_radii2 = sum_radii * sum_radii;
            if (dist2 < sum_radii2) {
              is_allowed = false;
              break;
            }
          }
        }
        if (is_allowed) {
          rescue_allowed[nb_allowed].dx = ix;
          rescue_allowed[nb_allowed].dy = iy;
          nb_allowed++;
        }
      }
    }
  } else {
    for (int ix = -search_zone_rescue.left; ix <= search_zone_rescue.right; ix++) {
      for (int iy = -search_zone_rescue.up; iy <= search_zone_rescue.down; iy++) {
        rescue_allowed[nb_allowed].dx = ix;
        rescue_allowed[nb_allowed].dy = iy;
        nb_allowed++;
      }
    }
  }

  for (drot = -rotmax; drot <= rotmax; drot += search_zone_rescue.inc_rot) {
    total_rot = refrot + grain[igrain].drot + drot;
    c = cos(total_rot);
    s = sin(total_rot);
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      rotate_pixel_pattern(igrain, i, c, s, &xpixel1, &ypixel1);
      grain[igrain].pattern1_rotated[i].dx = xpixel1;
      grain[igrain].pattern1_rotated[i].dy = ypixel1;
    }

    for (i_allowed = 0; i_allowed < nb_allowed; i_allowed++) {
      upix = rescue_allowed[i_allowed].dx;
      vpix = rescue_allowed[i_allowed].dy;
      mean1 = 0.0;
      C0C1 = C1C1 = 0.0;
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
        ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
        mean1 += image[im_index_current][xpixel1][ypixel1];
      }
      mean1 /= (double)grain[igrain].pattern.size();
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel0 = grain[igrain].pattern0_rotated[i].dx + refcoordx;
        ypixel0 = grain[igrain].pattern0_rotated[i].dy + refcoordy;
        xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
        ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
        C0C1 += (image[im_index_ref][xpixel0][ypixel0] - grain[igrain].mean0) *
                (image[im_index_current][xpixel1][ypixel1] - mean1);
        C1C1 +=
            (image[im_index_current][xpixel1][ypixel1] - mean1) * (image[im_index_current][xpixel1][ypixel1] - mean1);
      }

      NCC_test = C0C1 / sqrt(grain[igrain].C0C0 * C1C1);

      if (best_NCC < NCC_test) {
        best_NCC = NCC_test;
        best_upix = (double)upix;
        best_vpix = (double)vpix;
        best_drot = drot;
      }

    } // Loop i_allowed
  }   // Loop rot

  grain[igrain].NCC_rescue = best_NCC;

  if (best_NCC < NCC_min_super) {
#pragma omp critical
    { to_be_super_rescued[num_to_be_super_rescued++] = igrain; }
  } else {
    grain[igrain].upix = best_upix - rest_dx;
    grain[igrain].vpix = best_vpix - rest_dy;
    grain[igrain].rot_inc = best_drot;
    grain[igrain].dx = grain_dx + best_upix;
    grain[igrain].dy = grain_dy + best_vpix;
    grain[igrain].drot += best_drot;
  }
}

void follow_pattern_super_rescue_pixel(int igrain) {
  int refcoordx = grain[igrain].refcoord_xpix;
  int refcoordy = grain[igrain].refcoord_ypix;
  double refrot = grain[igrain].refrot;
  double c, s;
  int xpixel0, ypixel0;

  int grain_dx = nearest(grain[igrain].dx);
  int grain_dy = nearest(grain[igrain].dy);
  double rest_dx = grain[igrain].dx - grain_dx;
  double rest_dy = grain[igrain].dy - grain_dy;

  double drot;    // Increment test pour l'angle de rotation
  int upix, vpix; // Increment test pour la translation
  double mean1, C0C1, C1C1;

  int i_allowed;
  bool is_allowed;

  int xpixel1, ypixel1;
  double total_rot;

  double NCC_test, best_NCC = 0.0;
  double best_drot = 0.0;
  int best_upix = 0.0, best_vpix = 0.0;
  double rotmax = fabs(search_zone_super_rescue.inc_rot * search_zone_super_rescue.num_rot);

  // We define a list of relative coordinates used to search igrain in the image 1.
  // To do that, the list of neighbour of igrain is used
  const int nbmax = (search_zone_super_rescue.right + search_zone_super_rescue.left + 1) *
                    (search_zone_super_rescue.up + search_zone_super_rescue.down + 1);
  vector<relative_coord_type> super_rescue_allowed(nbmax);

  int radi_i = (int)grain[igrain].radius_pix;

  int nb_allowed = 0;

  if (use_neighbour_list) {
    for (int ix = -search_zone_super_rescue.left; ix <= search_zone_super_rescue.right; ix++) {
      for (int iy = -search_zone_super_rescue.up; iy <= search_zone_super_rescue.down; iy++) {
        is_allowed = true;
        int xcentre_i = refcoordx + grain_dx + ix;
        int ycentre_i = refcoordy + grain_dy + iy;
        for (int i_neigh = 0; i_neigh < grain[igrain].num_neighbour; i_neigh++) {
          int jgrain = grain[igrain].neighbour[i_neigh];
          double jgrain_NCC = max(grain[jgrain].NCC, grain[jgrain].NCC_rescue);
          if (jgrain_NCC > NCC_min_super) {
            int xcentre_j = grain[jgrain].refcoord_xpix + nearest(grain[jgrain].dx);
            int ycentre_j = grain[jgrain].refcoord_ypix + nearest(grain[jgrain].dy);
            int radi_j = (int)grain[jgrain].radius_pix;
            double dstx = xcentre_j - xcentre_i;
            double dsty = ycentre_j - ycentre_i;
            double dist2 = dstx * dstx + dsty * dsty;
            double sum_radii = (radi_i + radi_j) * 0.9; // only 90% of the radius declared
            double sum_radii2 = sum_radii * sum_radii;
            if (dist2 < sum_radii2) {
              is_allowed = false;
              break;
            }
          }
        }
        if (is_allowed) {
          super_rescue_allowed[nb_allowed].dx = ix;
          super_rescue_allowed[nb_allowed].dy = iy;
          nb_allowed++;
        }
      }
    }
  } else {
    for (int ix = -search_zone_super_rescue.left; ix <= search_zone_super_rescue.right; ix++) {
      for (int iy = -search_zone_super_rescue.up; iy <= search_zone_super_rescue.down; iy++) {
        super_rescue_allowed[nb_allowed].dx = ix;
        super_rescue_allowed[nb_allowed].dy = iy;
        nb_allowed++;
      }
    }
  }

  for (drot = -rotmax; drot <= rotmax; drot += search_zone_super_rescue.inc_rot) {
    total_rot = refrot + grain[igrain].drot + drot;
    c = cos(total_rot);
    s = sin(total_rot);
    for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
      rotate_pixel_pattern(igrain, i, c, s, &xpixel1, &ypixel1);
      grain[igrain].pattern1_rotated[i].dx = xpixel1;
      grain[igrain].pattern1_rotated[i].dy = ypixel1;
    }

    for (i_allowed = 0; i_allowed < nb_allowed; i_allowed++) {
      upix = super_rescue_allowed[i_allowed].dx;
      vpix = super_rescue_allowed[i_allowed].dy;
      mean1 = 0.0;
      C0C1 = C1C1 = 0.0;
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
        ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
        mean1 += (double)image[im_index_current][xpixel1][ypixel1];
      }
      mean1 /= (double)grain[igrain].pattern.size();
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel0 = grain[igrain].pattern0_rotated[i].dx + refcoordx;
        ypixel0 = grain[igrain].pattern0_rotated[i].dy + refcoordy;
        xpixel1 = grain[igrain].pattern1_rotated[i].dx + refcoordx + upix + grain_dx;
        ypixel1 = grain[igrain].pattern1_rotated[i].dy + refcoordy + vpix + grain_dy;
        C0C1 += ((double)image[im_index_ref][xpixel0][ypixel0] - grain[igrain].mean0) *
                ((double)image[im_index_current][xpixel1][ypixel1] - mean1);
        C1C1 += ((double)image[im_index_current][xpixel1][ypixel1] - mean1) *
                ((double)image[im_index_current][xpixel1][ypixel1] - mean1);
      }

      NCC_test = C0C1 / sqrt(grain[igrain].C0C0 * C1C1);

      if (best_NCC < NCC_test) {
        best_NCC = NCC_test;
        best_upix = (double)upix;
        best_vpix = (double)vpix;
        best_drot = drot;
      }
    }

  } // Loop for drot

  grain[igrain].NCC_rescue = best_NCC;

  if (best_NCC > NCC_min_super) { // 5 Juillet 2016. Si le super rescue est sans effet, alors on ne change pas
                                  // l'increment de deplacement
    grain[igrain].upix = best_upix - rest_dx; // (avril 2017) et les increments aussi...
    grain[igrain].vpix = best_vpix - rest_dy;
    grain[igrain].rot_inc = best_drot;

    grain[igrain].dx = grain_dx + best_upix;
    grain[igrain].dy = grain_dy + best_vpix;
    grain[igrain].drot += best_drot;
  }
}

double NCC_to_minimize_xyR(vector<double> &X) {

#if defined(_OPENMP)
  int thread_id = omp_get_thread_num();
#else
  int thread_id = 0;
#endif

  int igrain = igrain_of_thread[thread_id];

  double inv_num_pix = 1.0 / (double)(grain[igrain].pattern.size());
  int xpix0 = grain[igrain].refcoord_xpix;
  int ypix0 = grain[igrain].refcoord_ypix;
  double mean0 = grain[igrain].mean0;
  double C0C0 = grain[igrain].C0C0;

  double xpix1 = (double)grain[igrain].refcoord_xpix + (double)grain[igrain].dx + X[0];
  double ypix1 = (double)grain[igrain].refcoord_ypix + (double)grain[igrain].dy + X[1];

  double rot = grain[igrain].drot + X[2]; // from refrot which can be non-zero
  double c1 = cos(rot), s1 = sin(rot);

  vector<double> interpol_values;

  // Compute mean values
  double xsubpixel1, ysubpixel1;
  double mean1 = 0.0, int_grey_val;
  double x, y;

  for (size_t i = 0; i < grain[igrain].pattern0_rotated.size(); i++) {
    x = grain[igrain].pattern0_rotated[i].dx;
    y = grain[igrain].pattern0_rotated[i].dy;

    // Attention y vers le bas, le sens trigo doit etre inverse
    xsubpixel1 = xpix1 + c1 * (double)x + s1 * (double)y;
    ysubpixel1 = ypix1 - s1 * (double)x + c1 * (double)y;

    interpol_values.push_back(int_grey_val = IMAGE_INTERPOLATOR->getValue(im_index_current, xsubpixel1, ysubpixel1));
    mean1 += int_grey_val;
  }
  mean1 *= inv_num_pix;

  // NCC computation
  double C0C1 = 0.0, C1C1 = 0.0;
  double diffC0, diffC1;
  int xpixel0, ypixel0;
  int num_interpol_values = 0;
  for (size_t i = 0; i < grain[igrain].pattern0_rotated.size(); i++) {
    x = grain[igrain].pattern0_rotated[i].dx;
    y = grain[igrain].pattern0_rotated[i].dy;

    diffC1 = interpol_values[num_interpol_values++] - mean1;
    xpixel0 = xpix0 + x;
    ypixel0 = ypix0 + y;

    diffC0 = (double)image[im_index_ref][xpixel0][ypixel0] - mean0;
    C0C1 += diffC0 * diffC1;
    C1C1 += diffC1 * diffC1;
  }

  double NCC = C0C1 / sqrt(C0C0 * C1C1);

  double W = 1.0;
  double subpix_displacement = sqrt(X[0] * X[0] + X[1] * X[1]);
  if (subpix_displacement > subpix_displacement_max) {
    W = overflow_stiffness * (subpix_displacement - subpix_displacement_max) + 1.0;
  }

  return (1.0 - NCC) * W;
}

void follow_pattern_subpixel_xyR(int igrain) {
  vector<double> X(3); // Vector that hold the parameters to be optimized
  vector<double> DX(3);

  // ___ Initial guess of the subpixel part of displacement and rotation
  X[0] = X[1] = X[2] = 0.0;

  // ___ Initial direction (to search minimum in powell method)
  DX[0] = initial_direction_dx;
  DX[1] = initial_direction_dy;
  DX[2] = initial_direction_dr;

  // ___ Preparation of the minimization
#if defined(_OPENMP)
  int thread_id = omp_get_thread_num();
#else
  int thread_id = 0;
#endif

  igrain_of_thread[thread_id] = igrain;

  // ___ Minimization
  Powell<double(vector<double> &)> powell(NCC_to_minimize_xyR, subpix_tol);
  X = powell.minimize(X, DX);

  // ___ Update
  grain[igrain].NCC_subpix = 1.0 - powell.fret;

  grain[igrain].upix += X[0];
  grain[igrain].vpix += X[1];
  grain[igrain].rot_inc += X[2];

  grain[igrain].dx += X[0];
  grain[igrain].dy += X[1];
  grain[igrain].drot += X[2];
}

void find_neighbours(int igrain) {
  double x_igrain, y_igrain;
  double x_jgrain, y_jgrain;
  double squared_dist_pix = neighbour_dist_pix * neighbour_dist_pix;
  double dist2;

  x_igrain = grain[igrain].refcoord_xpix + grain[igrain].dx;
  y_igrain = grain[igrain].refcoord_ypix + grain[igrain].dy;

  grain[igrain].num_neighbour = 0;
  for (int jgrain = 0; jgrain < num_grains; jgrain++) {
    x_jgrain = grain[jgrain].refcoord_xpix + grain[jgrain].dx;
    y_jgrain = grain[jgrain].refcoord_ypix + grain[jgrain].dy;
    dist2 = (x_jgrain - x_igrain) * (x_jgrain - x_igrain) + (y_jgrain - y_igrain) * (y_jgrain - y_igrain);
    if (jgrain != igrain && dist2 <= squared_dist_pix) {
      grain[igrain].neighbour.push_back(jgrain);
    }
  }
}

int init(int argc, char *argv[]) {

  if (argc > 2) {
    fprintf(stderr, "Usage: %s command_file\n", argv[0]);
    fprintf(stderr, "Type %s -h for help\n", argv[0]);
    exit(0);
  }

  if (argc == 2) {
    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
      header();
      cout << endl;
      cout << "Help can be found in the folder 'examples' or 'doc'" << endl;
      cout << "  -v, --version   Display information about the soft" << endl;
      cout << "  -h, --help      This help" << endl;
      cout << "  -s, --sizes     Print size of common types (in bits) for this computer" << endl;
      cout << "  -g, --generate  Generate synthetic images for accuracy tests" << endl;
      cout << "  -d, --dialog    Use tracker to interactively extract information" << endl;
      cout << endl;
      exit(0);
    } else if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
      header();
      exit(0);
    } else if (strcmp(argv[1], "-s") == 0 || strcmp(argv[1], "--sizes") == 0) {
      cout << endl;
      cout << "unsigned char....." << sizeof(unsigned char) * 8 << " bits" << endl;
      cout << "unsigned short...." << sizeof(unsigned short) * 8 << " bits" << endl;
      cout << "unsigned int......" << sizeof(unsigned int) * 8 << " bits" << endl;
      cout << "float............." << sizeof(float) * 8 << " bits" << endl;
      cout << "double............" << sizeof(double) * 8 << " bits" << endl;
      cout << "uint16_t.........." << sizeof(uint16_t) * 8 << " bits" << endl;
      cout << "uint16_t..........." << sizeof(uint16_t) * 8 << " bits" << endl;
      cout << endl;
      exit(0);
    } else if (strcmp(argv[1], "-g") == 0 || strcmp(argv[1], "--generate") == 0) {
      int number_images = 0;
      double defx = 0.0, defy = 0.0, transx_pix = 0.0, transy_pix = 0.0, rot_deg = 0.0;
      double defx_inc = 0.0, defy_inc = 0.0, transx_pix_inc = 0.0, transy_pix_inc = 0.0, rot_deg_inc = 0.0;
      int w, h, octaves;
      double zoom, p;

      cout << endl;
      cout << "___Synthese of images." << endl;
      cout << "Number of images generated               : ";
      cin >> number_images;
      cout << "width (integer)                          : ";
      cin >> w;
      cout << "heigh (integer)                          : ";
      cin >> h;
      cout << endl;
      cout << "___Perlin noise. (octaves = 6, persistence = 0.5 and zoom = 5.0 seem to be ok)" << endl;
      cout << "octaves (integer)                        : ";
      cin >> octaves;
      cout << "persistence (double)                     : ";
      cin >> p;
      cout << "zoom (double)                            : ";
      cin >> zoom;
      cout << "___Transformations. (positive values correspond to left->right for x, top->bottom for y, "
              "counterclockwise for rotation)"
           << endl;
      cout << "x-deformation increment (-, double) : ";
      cin >> defx_inc;
      cout << "y-deformation increment (-, double) : ";
      cin >> defy_inc;
      cout << "x-translation increment (pixels, double) : ";
      cin >> transx_pix_inc;
      cout << "y-translation increment (pixels, double) : ";
      cin >> transy_pix_inc;
      cout << "rotation increment (degrees, double)     : ";
      cin >> rot_deg_inc;

      // LOG in a file
      ofstream file("synth_log.txt");
      file << "___Synthese of images." << endl;
      file << "Number of images generated               : " << number_images << endl;
      file << "width (integer)                          : " << w << endl;
      file << "heigh (integer)                          : " << h << endl;
      file << endl;
      file << "___Perlin noise." << endl;
      file << "octaves (integer)                        : " << octaves << endl;
      file << "persistence (double)                     : " << p << endl;
      file << "zoom (double)                            : " << zoom << endl;
      file << "___Transformations." << endl;
      file << "x-deformation increment (pixels, double) : " << defx_inc << endl;
      file << "y-deformation increment (pixels, double) : " << defy_inc << endl;
      file << "x-translation increment (pixels, double) : " << transx_pix_inc << endl;
      file << "y-translation increment (pixels, double) : " << transy_pix_inc << endl;
      file << "rotation increment (degrees, double)     : " << rot_deg_inc << endl;

      // image 0 non transformée
      generate_synthetic_images(w, h, zoom, octaves, p, defx, defy, transx_pix, transy_pix, rot_deg);
      for (int im = 1; im <= number_images; im++) {
        defx += defx_inc;
        defy += defy_inc;
        transx_pix += transx_pix_inc;
        transy_pix += transy_pix_inc;
        rot_deg += rot_deg_inc;

        generate_synthetic_images(w, h, zoom, octaves, p, defx, defy, transx_pix, transy_pix, rot_deg);
      }
      exit(0);
    } else if (strcmp(argv[1], "-d") == 0 || strcmp(argv[1], "--dialog") == 0) {
      dialog();
    } else {
      header();
      read_data(argv[1]);
    }
  } else { // tracker has been invocked by double-clics or without arguments
    header();
    read_data("commands.txt");
  }

  g_logfile.open("tracking.log");

  if (rotations == 0) { // to be replaced by if (search_zone_super_rescue.inc_rot == 0.0 || ...) (FIXME)
    LOG("Rotation are NOT TRACKED!" << endl);
    // Il faut choisir une valeur non nulle pour inc_rot sinon -> boucle infinie
    search_zone_super_rescue.inc_rot = 1.0;
    search_zone_rescue.inc_rot = 1.0;
    search_zone.inc_rot = 1.0;

    search_zone_super_rescue.num_rot = 0;
    search_zone_rescue.num_rot = 0;
    search_zone.num_rot = 0;
  }

#ifdef _OPENMP
  omp_set_num_threads(wanted_num_threads);
#endif
  // fprintf(stdout, "Number of threads: %d\n", wanted_num_threads);
  LOG("Number of threads: " << wanted_num_threads << endl);

  return 1;
}

int read_gmsh(const char *name) {
  ifstream meshFile(name);
  if (!meshFile) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  string token;
  char not_read[150];
  meshFile >> token;

  while (meshFile) {
    if (token == "$Nodes") {
      meshFile >> num_grains;
      SHOW(num_grains);

      if (!grain.empty())
        grain.clear();
      grain.resize(num_grains);
      if (!to_be_rescued.empty())
        to_be_rescued.clear();
      to_be_rescued.resize(num_grains);
      if (!to_be_super_rescued.empty())
        to_be_super_rescued.clear();
      to_be_super_rescued.resize(num_grains);

      size_t num_node;
      double radius_pix_real = 10.0;
      double xpix_real, ypix_real, trash;
      for (int i = 0; i < num_grains; i++) {
        meshFile >> num_node >> xpix_real >> ypix_real >> trash;
        grain[i].refrot = 0.0;
        grain[i].radius_pix = radius_pix_real;
        grain[i].refcoord_xpix = (int)floor(xpix_real);
        grain[i].dx = grain[i].upix = xpix_real - floor(xpix_real);
        grain[i].refcoord_ypix = (int)floor(ypix_real);
        grain[i].dy = grain[i].vpix = ypix_real - floor(ypix_real);
      }
    }

    if (token == "$Elements") {
      mesh.clear();
      size_t nbElements;
      size_t num_element, element_type, nbTags;
      triangle T;
      size_t t = 0;

      meshFile >> nbElements;
      for (size_t e = 0; e < nbElements; ++e) {
        meshFile >> num_element >> element_type;
        // triangle
        if (element_type != 2) {
          meshFile.getline(not_read, 150);
          continue;
        }

        meshFile >> nbTags;
        // the third tag is the number of a mesh partition to which the element belongs
        size_t tag;
        for (size_t tg = 0; tg < nbTags; ++(tg)) {
          meshFile >> tag;
        }

        meshFile >> T.i0 >> T.i1 >> T.i2;

        // nodeId has 0-offset
        // (0 in C/C++ corresponds to 1 in the file)
        T.i0 -= 1;
        T.i1 -= 1;
        T.i2 -= 1;

        mesh.push_back(T);
        ++t;
      }
    }

    if (token == "$EndElements")
      break;

    meshFile >> token;
  }
  return 1;
}

// Lecture des positions des grains "tracked"
int read_grains(const char *name, bool restart) {
  ifstream grain_file(name);
  if (!grain_file) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  grain_file >> num_grains;
  SHOW(num_grains);

  if (!grain.empty())
    grain.clear();
  grain.resize(num_grains);
  if (!to_be_rescued.empty())
    to_be_rescued.clear();
  to_be_rescued.resize(num_grains);
  if (!to_be_super_rescued.empty())
    to_be_super_rescued.clear();
  to_be_super_rescued.resize(num_grains);

  if (restart) { // File with the same format as 'dic_out_x.txt'
    for (int i = 0; i < num_grains; i++) {
      grain_file >> grain[i].refcoord_xpix >> grain[i].refcoord_ypix >> grain[i].refrot >> grain[i].radius_pix >>
          grain[i].dx >> grain[i].dy >> grain[i].drot >> grain[i].upix >> grain[i].vpix >> grain[i].rot_inc >>
          grain[i].NCC >> grain[i].NCC_rescue >> grain[i].NCC_subpix;
    }
  } else { // Format plus simple
    double radius_pix_real;
    double xpix_real, ypix_real;
    for (int i = 0; i < num_grains; i++) {
      grain_file >> xpix_real // 1
          >> ypix_real        // 2
          >> grain[i].refrot  // 3
          >> radius_pix_real; // 4
      grain[i].radius_pix = radius_pix_real;

      // Attention, si les coordonnées sont des reels, grain.dx est alors non nul
      // et cela ajoute un déplacement initial artificiel
      grain[i].refcoord_xpix = (int)floor(xpix_real);
      grain[i].dx = grain[i].upix = xpix_real - floor(xpix_real);
      grain[i].refcoord_ypix = (int)floor(ypix_real);
      grain[i].dy = grain[i].vpix = ypix_real - floor(ypix_real);
    }
  }

  return 1;
}

void make_grid(double xmin, double xmax, double ymin, double ymax, int nx, int ny, int aleaMax) {
  int dx = (int)floor((xmax - xmin) / (double)(nx - 1));
  int dy = (int)floor((ymax - ymin) / (double)(ny - 1));

  // Allocation memoire pour grain et to_be_rescued
  num_grains = nx * ny;
  if (!grain.empty())
    grain.clear();
  grain.resize(num_grains);
  if (!to_be_rescued.empty())
    to_be_rescued.clear();
  to_be_rescued.resize(num_grains);
  if (!to_be_super_rescued.empty())
    to_be_super_rescued.clear();
  to_be_super_rescued.resize(num_grains);

  cout << "Create grid with " << num_grains << " points" << endl;

  int i = 0;
  for (int iy = 0; iy < ny; iy++) {
    for (int ix = 0; ix < nx; ix++) {
      double angle = (double)rand() / (double)RAND_MAX * (2.0 * M_PI); // un angle entre 0 et 2 pi
      double alea = (double)rand() / (double)RAND_MAX * aleaMax;
      grain[i].refcoord_xpix = xmin + ix * dx + (int)nearest(alea * cos(angle));
      grain[i].refcoord_ypix = ymin + iy * dy + (int)nearest(alea * sin(angle));
      grain[i].refrot = 0.0;
      grain[i].radius_pix = 2.0;
      i++;
    }
  }
}

void make_grid_circular(double x_center, double y_center, double radius, double angle_inc) {
  int i = 0;
  double Pi = 4.0 * atan(1.0);

  // Allocation memoire pour grain et to_be_rescued
  num_grains = (int)(2.0 * Pi / angle_inc + 1.0);
  if (!grain.empty())
    grain.clear();
  grain.resize(num_grains);
  if (!to_be_rescued.empty())
    to_be_rescued.clear();
  to_be_rescued.resize(num_grains);
  if (!to_be_super_rescued.empty())
    to_be_super_rescued.clear();
  to_be_super_rescued.resize(num_grains);

  for (double angle = 0; angle < 2.0 * Pi; angle += angle_inc) {
    double x = x_center + radius * cos(angle);
    double y = y_center + radius * sin(angle);
    grain[i].refcoord_xpix = nearest(x);
    grain[i].refcoord_ypix = nearest(y);
    grain[i].refrot = 0.0;
    grain[i].radius_pix = 1.0;
    i++;
  }
}

int save_grains(const char *name, int num, bool simpleVersion) {
  ofstream grain_file_out(name);
  if (!grain_file_out) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  if (simpleVersion) {
    grain_file_out << grain.size() << endl;
    for (size_t i = 0; i < grain.size(); i++) {
      grain_file_out << grain[i].refcoord_xpix + grain[i].dx << ' ' << grain[i].refcoord_ypix + grain[i].dy << ' '
                     << grain[i].refrot << ' ' << grain[i].radius_pix << '\n';
    }
  } else {
    grain_file_out << num_grains << endl;
    for (int i = 0; i < num_grains; i++) {
      grain_file_out << scientific << setprecision(15) << setw(6) << left << grain[i].refcoord_xpix << ' ' // 1
                     << setw(6) << left << grain[i].refcoord_ypix << ' '                                   // 2
                     << setw(23) << right << grain[i].refrot << ' '                                        // 3
                     << setw(23) << right << grain[i].radius_pix << ' '                                    // 4
                     << setw(23) << right << grain[i].dx << ' '                                            // 5
                     << setw(23) << right << grain[i].dy << ' '                                            // 6
                     << setw(23) << right << grain[i].drot << ' '                                          // 7
                     << setw(23) << right << grain[i].upix << ' '                                          // 8
                     << setw(23) << right << grain[i].vpix << ' '                                          // 9
                     << setw(23) << right << grain[i].rot_inc << ' '                                       // 10
                     << setw(23) << right << grain[i].NCC << ' '                                           // 11
                     << setw(23) << right << grain[i].NCC_rescue << ' '                                    // 12
                     << setw(23) << right << grain[i].NCC_subpix << '\n';                                  // 13
    }

    grain_file_out << endl;

    char fname[256];
    sprintf(fname, image_name, iref);
    grain_file_out << "# from image: " << fname << " --dateTime: " << imageData[im_index_ref].dateTime ; //<< "\n";
    sprintf(fname, image_name, num);
    grain_file_out << "#   to image: " << fname << " --dateTime: " << imageData[im_index_current].dateTime ; //<< "\n";
    //grain_file_out << fname << endl;
  }

  return 1;
}

int save_grains(int num) {
  char name[256];
  sprintf(name, "dic_out_%d.txt", num);

  return save_grains(name, num);
}

void make_rect_pattern(int igrain, int half_width_pix, int half_height_pix) {
  // Reserver memoire pour pattern, pattern0_rotation et pattern1_rotation
  int dim = (2 * half_width_pix + 1) * (2 * half_height_pix + 1);
  grain[igrain].pattern.resize(dim);
  grain[igrain].pattern0_rotated.resize(dim);
  grain[igrain].pattern1_rotated.resize(dim);

  int i = 0;
  for (int dy = -half_height_pix; dy <= half_height_pix; ++dy) {
    for (int dx = -half_width_pix; dx <= half_width_pix; ++dx) {
      grain[igrain].pattern[i].dx = dx;
      grain[igrain].pattern[i].dy = dy;
      i++;
    }
  }
}

void make_circ_pattern(int igrain, int radius_pix) {
  // Reserver memoire pour pattern, pattern0_rotation et pattern1_rotation
  int dim = 4 * radius_pix * radius_pix; // surdim !
  grain[igrain].pattern.resize(dim);
  grain[igrain].pattern0_rotated.resize(dim);
  grain[igrain].pattern1_rotated.resize(dim);
  int i = 0;
  double dst;
  for (int dy = -radius_pix; dy <= radius_pix; ++dy) {
    for (int dx = -radius_pix; dx <= radius_pix; ++dx) {
      dst = sqrt((double)dx * (double)dx + (double)dy * (double)dy);
      if (dst <= (double)radius_pix) {
        grain[igrain].pattern[i].dx = dx;
        grain[igrain].pattern[i].dy = dy;
        i++;
      }
    }
  }

  grain[igrain].pattern.resize(i);
  grain[igrain].pattern0_rotated.resize(i);
  grain[igrain].pattern1_rotated.resize(i);
}

void make_ring_pattern(int igrain, int radius_IN_pix, int radius_OUT_pix) {
  // Reserver memoire pour pattern, pattern0_rotation et pattern1_rotation
  int dim = 4 * radius_OUT_pix * radius_OUT_pix; // surdim !
  grain[igrain].pattern.resize(dim);
  grain[igrain].pattern0_rotated.resize(dim);
  grain[igrain].pattern1_rotated.resize(dim);

  int i = 0;
  double dst;
  for (int dy = -radius_OUT_pix; dy <= radius_OUT_pix; ++dy) {
    for (int dx = -radius_OUT_pix; dx <= radius_OUT_pix; ++dx) {
      dst = sqrt((double)dx * (double)dx + (double)dy * (double)dy);
      if (dst <= (double)radius_OUT_pix && dst >= (double)radius_IN_pix) {
        grain[igrain].pattern[i].dx = dx;
        grain[igrain].pattern[i].dy = dy;
        i++;
      }
    }
  }

  grain[igrain].pattern.resize(i);
  grain[igrain].pattern0_rotated.resize(i);
  grain[igrain].pattern1_rotated.resize(i);
}

// Format:
// number_of_patterns
// for each pattern {
// 	number_of_relative_coords
// 	for each coord {
//		xrel yrel
// 	}
// }
void make_custom_pattern(const char *name) {
  ifstream pattern_file(name);
  if (!pattern_file) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  int num_patterns;
  pattern_file >> num_patterns;
  if (num_patterns != num_grains) {
    fprintf(stderr, "The number of patterns is not equal to the number of grains!\n");
    fprintf(stderr, "num_patterns : %d\n", num_patterns);
    fprintf(stderr, "num_grains : %d\n", num_grains);
    exit(0);
  }
  fprintf(stdout, "Number of patterns %d\n", num_patterns);

  int num_coords = 0;
  int dx, dy;
  for (int igrain = 0; igrain < num_patterns; igrain++) {
    pattern_file >> num_coords;
    // Reserver memoire pour pattern, pattern0_rotation et pattern1_rotation
    grain[igrain].pattern.resize(num_coords);
    grain[igrain].pattern0_rotated.resize(num_coords);
    grain[igrain].pattern1_rotated.resize(num_coords);

    // grain[igrain].num_point_pattern = num_coords;
    for (int i = 0; i < num_coords; i++) {
      pattern_file >> dx >> dy;
      grain[igrain].pattern[i].dx = dx;
      grain[igrain].pattern[i].dy = dy;
    }
  }
}

void set_solidTranslation(double u, double v) {
  for (int igrain = 0; igrain < num_grains; igrain++) {
    grain[igrain].dx = grain[igrain].upix = u;
    grain[igrain].dy = grain[igrain].vpix = v;
  }
}

void mask_rect(int xmin, int xmax, int ymin, int ymax) {
  for (int igrain = 0; igrain < num_grains; igrain++) {
    if (grain[igrain].refcoord_xpix >= xmin && grain[igrain].refcoord_xpix <= xmax &&
        grain[igrain].refcoord_ypix >= ymin && grain[igrain].refcoord_ypix <= ymax)
      grain[igrain].masked = true;
  }
}

// This function reads the command file
int read_data(const char *name) {
  ifstream command_file(name);
  if (!command_file) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  string token;
  command_file >> token;

  while (command_file) {

    if (token[0] == '!')
      getline(command_file, token);
    else if (token == "image_name") {
      command_file >> image_name;
    } else if (token == "dic_name") {
      command_file >> dic_name;
    } else if (token == "RawImages") {
      command_file >> RawImages;
    } else if (token == "DemosaicModel") {
      command_file >> DemosaicModel;
    } else if (token == "rescaleGrayLevels") {
      command_file >> rescaleGrayLevels;
    } else if (token == "grains_to_follow" || token == "start_file") {
      char filename[256];
      command_file >> filename;
      // grain_positions_command = "grains_to_follow " + filename;
      read_grains(filename, false);
    } else if (token == "gmsh_file") {
      char filename[256];
      command_file >> filename;
      read_gmsh(filename);
    } else if (token == "restart_file") {
      char filename[256];
      command_file >> filename;
      // grain_positions_command = "restart_file " + filename;
      read_grains(filename, true);
    } else if (token == "reordering") {
      string order;
      command_file >> order;
      if (order == "increasingRadius") {
        sort(grain.begin(), grain.end(), increasingRadius_order());
        save_grains("increasingRadius.txt", 0, true);
        cout << "Grains have been reordered with increasing radius in file increasingRadius.txt" << endl;
        exit(0);
      } else if (order == "decreasingRadius") {
        sort(grain.begin(), grain.end(), decreasingRadius_order());
        save_grains("decreasingRadius.txt", 0, true);
        cout << "Grains have been reordered with decreasing radius in file decreasingRadius.txt" << endl;
        exit(0);
      } else if (order == "increasingHeight") {
        sort(grain.begin(), grain.end(), increasingHeight_order());
        save_grains("increasingHeight.txt", 0, true);
        cout << "Grains have been reordered with increasing height in file increasingHeight.txt" << endl;
        exit(0);
      } else if (order == "mixed") {
        std::random_shuffle(grain.begin(), grain.end());
        save_grains("mixed.txt", 0, true);
        cout << "Grains have been mixed in file mixed.txt" << endl;
        exit(0);
      }
    } else if (token == "make_grid") {
      double xmin, xmax, ymin, ymax;
      int nx, ny;
      int aleaMax;
      command_file >> xmin >> xmax >> ymin >> ymax >> nx >> ny >> aleaMax;
      // grain_positions_command = "make_grid " + string(xmin) + " ..." ;
      make_grid(xmin, xmax, ymin, ymax, nx, ny, aleaMax);
    } else if (token == "make_grid_circular") {
      double x_center, y_center, radius, angle_inc;
      command_file >> x_center >> y_center >> radius >> angle_inc;
      // grain_positions_command = "make_grid_circular " + string(x_center) + " ..." ;
      make_grid_circular(x_center, y_center, radius, angle_inc);
    } else if (token == "iref") {
      command_file >> iref;
    } else if (token == "ibeg") {
      command_file >> ibeg;
    } else if (token == "iend") {
      command_file >> iend;
    } else if (token == "iinc") {
      command_file >> iinc;
    } else if (token == "iraz") {
      command_file >> iraz;
    } else if (token == "idelta") {
      command_file >> idelta;
    }

    else if (token == "wanted_num_threads") {
      int n;
      command_file >> n;
      wanted_num_threads = fabs(n);
    } else if (token == "procedure") {
      command_file >> procedure;
    } else if (token == "verbose_level") {
      command_file >> verbose_level;
    } else if (token == "rescue_level") {
      command_file >> rescue_level;
    } else if (token == "subpixel") {
      command_file >> subpixel;
    } else if (token == "rotations") {
      command_file >> rotations;
    }

    else if (token == "num_neighbour_max") {
      command_file >> num_neighbour_max;
    } else if (token == "neighbour_dist_pix") {
      command_file >> neighbour_dist_pix;
    } else if (token == "period_rebuild_neighbour_list") {
      command_file >> period_rebuild_neighbour_list;
    }

    else if (token == "image_interpolator") {
      command_file >> token;
      if (token == "linear") {
        IMAGE_INTERPOLATOR = &IMAGE_INTERPOLATOR_LINEAR;
      } else if (token == "cubic") {
        IMAGE_INTERPOLATOR = &IMAGE_INTERPOLATOR_CUBIC;
      } else if (token == "quintic") {
        IMAGE_INTERPOLATOR = &IMAGE_INTERPOLATOR_QUINTIC;
      }
    } else if (token == "interp_mask_size") {
      cout << "WARNING!! interp_mask_size is deprecated since version 0.5.0\n";
      int a;
      command_file >> a;
    } else if (token == "subpix_tol") {
      command_file >> subpix_tol;
    } else if (token == "initial_direction_dx") {
      command_file >> initial_direction_dx;
    } else if (token == "initial_direction_dy") {
      command_file >> initial_direction_dy;
    } else if (token == "initial_direction_dr") {
      command_file >> initial_direction_dr;
    } else if (token == "overflow_stiffness") {
      command_file >> overflow_stiffness;
    } else if (token == "subpix_displacement_max") {
      command_file >> subpix_displacement_max;
    }

    else if (token == "make_images") {
      command_file >> make_images;
    } else if (token == "image_size") {
      command_file >> image_size;
      image_div = nearest(1. / image_size);
    } else if (token == "image_div") {
      command_file >> image_div;
      image_size = 1.0 / (double)image_div;
    } else if (token == "draw_angle") {
      command_file >> draw_angle;
    } else if (token == "draw_disp") {
      command_file >> draw_disp;
    } else if (token == "draw_rescued") {
      command_file >> draw_rescued;
    }

    else if (token == "search_zone.left") {
      command_file >> search_zone.left;
    } else if (token == "search_zone.right") {
      command_file >> search_zone.right;
    } else if (token == "search_zone.up") {
      command_file >> search_zone.up;
    } else if (token == "search_zone.down") {
      command_file >> search_zone.down;
    } else if (token == "search_zone.inc_rot") {
      command_file >> search_zone.inc_rot;
    } else if (token == "search_zone.num_rot") {
      command_file >> search_zone.num_rot;
    }

    else if (token == "NCC_min") {
      command_file >> NCC_min;
    } else if (token == "search_zone_rescue.left") {
      command_file >> search_zone_rescue.left;
    } else if (token == "search_zone_rescue.right") {
      command_file >> search_zone_rescue.right;
    } else if (token == "search_zone_rescue.up") {
      command_file >> search_zone_rescue.up;
    } else if (token == "search_zone_rescue.down") {
      command_file >> search_zone_rescue.down;
    } else if (token == "search_zone_rescue.inc_rot") {
      command_file >> search_zone_rescue.inc_rot;
    } else if (token == "search_zone_rescue.num_rot") {
      command_file >> search_zone_rescue.num_rot;
    }

    else if (token == "NCC_min_super") {
      command_file >> NCC_min_super;
    } else if (token == "search_zone_super_rescue.left") {
      command_file >> search_zone_super_rescue.left;
    } else if (token == "search_zone_super_rescue.right") {
      command_file >> search_zone_super_rescue.right;
    } else if (token == "search_zone_super_rescue.up") {
      command_file >> search_zone_super_rescue.up;
    } else if (token == "search_zone_super_rescue.down") {
      command_file >> search_zone_super_rescue.down;
    } else if (token == "search_zone_super_rescue.inc_rot") {
      command_file >> search_zone_super_rescue.inc_rot;
    } else if (token == "search_zone_super_rescue.num_rot") {
      command_file >> search_zone_super_rescue.num_rot;
    } else if (token == "use_neighbour_list") {
      command_file >> use_neighbour_list;
    }

    else if (token == "targetRadiusPattern") {
      command_file >> targetRadiusPattern;
    } else if (token == "pattern") {
      command_file >> token;
      if (token == "circ") {
        int radius;
        command_file >> radius;
        // Same pattern for each grain
        for (int i = 0; i < num_grains; i++)
          if (targetRadiusPattern < 0.0 || grain[i].radius_pix == targetRadiusPattern)
            make_circ_pattern(i, radius);
      } else if (token == "grain_diameter") {
        double radius_reduction, min_radius;
        command_file >> radius_reduction >> min_radius;
        // circular pattern for each grain with a radius dependant to grain radius
        for (int i = 0; i < num_grains; i++) {
          int radius = (int)(grain[i].radius_pix * radius_reduction);
          if (radius < min_radius) {
            radius = min_radius;
          }
          if (targetRadiusPattern < 0.0 || grain[i].radius_pix == targetRadiusPattern)
            make_circ_pattern(i, radius);
        }
      } else if (token == "ring") {
        int radius_IN, radius_OUT;
        command_file >> radius_IN >> radius_OUT;
        // Same pattern for each grain
        for (int i = 0; i < num_grains; i++)
          if (targetRadiusPattern < 0.0 || grain[i].radius_pix == targetRadiusPattern)
            make_ring_pattern(i, radius_IN, radius_OUT);
      } else if (token == "file") {
        char filename[256];
        command_file >> filename;
        make_custom_pattern(filename);
      } else if (token == "rect") {
        int half_width_pix, half_height_pix;
        command_file >> half_width_pix >> half_height_pix;
        // Same pattern for each grain
        for (int i = 0; i < num_grains; i++)
          if (targetRadiusPattern < 0.0 || grain[i].radius_pix == targetRadiusPattern)
            make_rect_pattern(i, half_width_pix, half_height_pix);
      }
    }

    else if (token == "subpix_center_threshold") {
      command_file >> subpix_center_threshold;
    } else if (token == "subpix_center_dr0") {
      command_file >> subpix_center_dr0;
    } else if (token == "subpix_center_xstep") {
      command_file >> subpix_center_xstep;
    } else if (token == "subpix_center_ystep") {
      command_file >> subpix_center_ystep;
    } else if (token == "subpix_center_rstep") {
      command_file >> subpix_center_rstep;
    }

    else if (token == "image_numbers_corrDisto") {
      size_t nb;
      command_file >> nb;
      image_numbers_corrDisto.resize(nb);
      for (size_t i = 0; i < nb; i++) {
        command_file >> image_numbers_corrDisto[i];
      }
    } else if (token == "imposedDisplApprox_corrDisto") {
      imposed_displ_x_pix.resize(image_numbers_corrDisto.size());
      imposed_displ_y_pix.resize(image_numbers_corrDisto.size());
      imposed_displ_x_pix[0] = imposed_displ_y_pix[0] = 0.0;
      for (size_t i = 1; i < image_numbers_corrDisto.size(); i++) {
        command_file >> imposed_displ_x_pix[i] >> imposed_displ_y_pix[i];
      }
    } else if (token == "setSolidSolutionFromRef") {
      double dx, dy, drot;
      command_file >> dx >> dy >> drot;
      setSolidSolutionFromRef(dx, dy, drot);
    } else if (token == "ShowSolidMotionError") {
      command_file >> ShowSolidMotionError;
    }

    else if (token == "grid_image_name") {
      command_file >> grid_image_name;
    } else if (token == "grid_dim") {
      command_file >> nx_grid_disto >> ny_grid_disto;
    }

    else if (token == "fake_undistor") {
      command_file >> fake_undistor;
    }

    else if (token == "equiProj_dist_range") {
      command_file >> equiProj_dist_min >> equiProj_dist_max;
    }

    else if (token == "xc_corrDistor") {
      command_file >> disto_parameters[0];
    } else if (token == "yc_corrDistor") {
      command_file >> disto_parameters[1];
    } else if (token == "K1_corrDistor") {
      command_file >> disto_parameters[2];
    } else if (token == "K2_corrDistor") {
      command_file >> disto_parameters[3];
    } else if (token == "K3_corrDistor") {
      command_file >> disto_parameters[4];
    } else if (token == "P1_corrDistor") {
      command_file >> disto_parameters[5];
    } else if (token == "P2_corrDistor") {
      command_file >> disto_parameters[6];
    } else if (token == "P3_corrDistor") {
      command_file >> disto_parameters[7];
    }

    else if (token == "dxc_corrDistor") {
      command_file >> disto_parameters_perturb[0];
    } else if (token == "dyc_corrDistor") {
      command_file >> disto_parameters_perturb[1];
    } else if (token == "dK1_corrDistor") {
      command_file >> disto_parameters_perturb[2];
    } else if (token == "dK2_corrDistor") {
      command_file >> disto_parameters_perturb[3];
    } else if (token == "dK3_corrDistor") {
      command_file >> disto_parameters_perturb[4];
    } else if (token == "dP1_corrDistor") {
      command_file >> disto_parameters_perturb[5];
    } else if (token == "dP2_corrDistor") {
      command_file >> disto_parameters_perturb[6];
    } else if (token == "dP3_corrDistor") {
      command_file >> disto_parameters_perturb[7];
    }

    else if (token == "undistor") {
      string filename;
      command_file >> filename;
      string filename_to = "undist_" + filename;
      undistor_image(filename.c_str(), filename_to.c_str());
    } else if (token == "halfPatternQual") {
      command_file >> halfPatternQual;
    }

    else if (token == "colorMin") {
      command_file >> colorMin;
    } else if (token == "colorMax") {
      command_file >> colorMax;
    } else if (token == "nx_grid_visu") {
      command_file >> nx_grid_visu;
    } else if (token == "smoothDegreeF") {
      command_file >> smoothDegreeF;
    } else if (token == "visu_autoscale") {
      command_file >> visu_autoscale;
    } else if (token == "visu_alpha") {
      command_file >> visu_alpha;
    } else if (token == "visu_output") {
      command_file >> visu_output;
    } else if (token == "visu_mode") {
      command_file >> visu_mode;
      cout << "visu_mode = " << visu_mode << endl;
    } else if (token == "visu_draw_in_ref") {
      command_file >> visu_draw_in_ref;
    }

    else if (token == "set_solidTranslation") {
      double u, v;
      command_file >> u >> v;
      set_solidTranslation(u, v);
    } else if (token == "mask_rect") {
      int xmin, xmax, ymin, ymax;
      command_file >> xmin >> xmax >> ymin >> ymax;
      mask_rect(xmin, xmax, ymin, ymax);
    }

    else {
      fprintf(stdout, "Unknown token: %s\n", token.c_str());
    }

    command_file >> token;
  }

  return 1;
}

// This function writes the command file
void write_data(const char *name) {
  ofstream command_file(name);
  if (!command_file) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  // string token;
  // command_file >> token;

  command_file << "image_name " << image_name << endl;
  command_file << "RawImages " << RawImages << endl;
  command_file << "DemosaicModel " << DemosaicModel << endl;

  command_file << grain_positions_command << endl; //// TODO  !!!!!!!!!!!!!!!!

  command_file << "iref " << iref << endl;
  command_file << "ibeg " << ibeg << endl;
  command_file << "iend " << iend << endl;
  command_file << "iinc " << iinc << endl;
  command_file << "iraz " << iraz << endl;

  command_file << "wanted_num_threads " << wanted_num_threads << endl;
  command_file << "procedure " << procedure << endl;
  command_file << "verbose_level " << verbose_level << endl;
  command_file << "rescue_level " << rescue_level << endl;
  command_file << "subpixel " << subpixel << endl;
  command_file << "rotations " << rotations << endl;

  command_file << "num_neighbour_max " << num_neighbour_max << endl;
  command_file << "neighbour_dist_pix " << neighbour_dist_pix << endl;
  command_file << "period_rebuild_neighbour_list " << period_rebuild_neighbour_list << endl;

  command_file << "image_interpolator "
               << "linear" << endl; ////// TODO !!!!!!!!!!!!!!
  command_file << "subpix_tol " << subpix_tol << endl;
  command_file << "initial_direction_dx " << initial_direction_dx << endl;
  command_file << "initial_direction_dy " << initial_direction_dy << endl;
  command_file << "initial_direction_dr " << initial_direction_dr << endl;

  command_file << "make_images " << make_images << endl;
  command_file << "image_size " << image_size << endl;
  command_file << "draw_angle " << draw_angle << endl;
  command_file << "draw_disp " << draw_disp << endl;
  command_file << "draw_rescued " << draw_rescued << endl;

  command_file << "search_zone.left " << search_zone.left << endl;
  command_file << "search_zone.right " << search_zone.right << endl;
  command_file << "search_zone.up " << search_zone.up << endl;
  command_file << "search_zone.down " << search_zone.down << endl;
  command_file << "search_zone.inc_rot " << search_zone.inc_rot << endl;
  command_file << "search_zone.num_rot " << search_zone.num_rot << endl;

  command_file << "NCC_min " << NCC_min << endl;
  command_file << "search_zone_rescue.left " << search_zone_rescue.left << endl;
  command_file << "search_zone_rescue.right " << search_zone_rescue.right << endl;
  command_file << "search_zone_rescue.up " << search_zone_rescue.up << endl;
  command_file << "search_zone_rescue.down " << search_zone_rescue.down << endl;
  command_file << "search_zone_rescue.inc_rot " << search_zone_rescue.inc_rot << endl;
  command_file << "search_zone_rescue.num_rot " << search_zone_rescue.num_rot << endl;

  command_file << "NCC_min_super " << NCC_min_super << endl;
  command_file << "search_zone_super_rescue.left " << search_zone_super_rescue.left << endl;
  command_file << "search_zone_super_rescue.right " << search_zone_super_rescue.right << endl;
  command_file << "search_zone_super_rescue.up " << search_zone_super_rescue.up << endl;
  command_file << "search_zone_super_rescue.down " << search_zone_super_rescue.down << endl;
  command_file << "search_zone_super_rescue.inc_rot " << search_zone_super_rescue.inc_rot << endl;
  command_file << "search_zone_super_rescue.num_rot " << search_zone_super_rescue.num_rot << endl;
  command_file << "use_neighbour_list " << use_neighbour_list << endl;

  command_file << pattern_command << endl; // TODO !!!!!!!!

  command_file << "subpix_center_threshold " << subpix_center_threshold << endl;
  command_file << "subpix_center_dr0 " << subpix_center_dr0 << endl;
  command_file << "subpix_center_rstep " << subpix_center_rstep << endl;

  // command_file << "image_numbers_corrDisto " << XXXX << endl;
  // command_file << "imposedDisplApprox_corrDisto " << XXXX << endl;

  command_file << "xc_corrDistor " << disto_parameters[0] << endl;
  command_file << "yc_corrDistor " << disto_parameters[1] << endl;
  command_file << "K1_corrDistor " << disto_parameters[2] << endl;
  command_file << "K2_corrDistor " << disto_parameters[3] << endl;
  command_file << "K3_corrDistor " << disto_parameters[4] << endl;
  command_file << "P1_corrDistor " << disto_parameters[5] << endl;
  command_file << "P2_corrDistor " << disto_parameters[6] << endl;
  command_file << "P3_corrDistor " << disto_parameters[7] << endl;

  command_file << "dxc_corrDistor " << disto_parameters_perturb[0] << endl;
  command_file << "dyc_corrDistor " << disto_parameters_perturb[1] << endl;
  command_file << "dK1_corrDistor " << disto_parameters_perturb[2] << endl;
  command_file << "dK2_corrDistor " << disto_parameters_perturb[3] << endl;
  command_file << "dK3_corrDistor " << disto_parameters_perturb[4] << endl;
  command_file << "dP1_corrDistor " << disto_parameters_perturb[5] << endl;
  command_file << "dP2_corrDistor " << disto_parameters_perturb[6] << endl;
  command_file << "dP3_corrDistor " << disto_parameters_perturb[7] << endl;

  command_file << "! END OF FILE";
}

// PROCEDURES:

#include "distortion_correction.cpp"
#include "grayLevelAnalysis.cpp"
#include "pattern_quality.cpp"
#include "post_process.cpp"
#include "subpixel_centers.cpp"
#include "synthetic_images.cpp"
#include "visu_process.cpp"
//#include "solid_motion_error.cpp"
