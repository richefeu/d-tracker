// ------------------------------------------
// D-TRACKER, December 2011-2026
// Laboratoire 3SR, Université Grenoble Alpes
// ------------------------------------------

#define TRACKER_VERSION "0.7.0"

#include "tracker.hpp"

void header() {
  // http://www.patorjk.com/software/taag
  // Lean
  std::cout << "\n";
  std::cout << "_/_/_/_/_/  _/_/_/      _/_/      _/_/_/  _/    _/  _/_/_/_/  _/_/_/ " << "\n";
  std::cout << "   _/      _/    _/  _/    _/  _/        _/  _/    _/        _/    _/" << "\n";
  std::cout << "  _/      _/_/_/    _/_/_/_/  _/        _/_/      _/_/_/    _/_/_/   " << "\n";
  std::cout << " _/      _/    _/  _/    _/  _/        _/  _/    _/        _/    _/  " << "\n";
  std::cout << "_/      _/    _/  _/    _/    _/_/_/  _/    _/  _/_/_/_/  _/    _/   " << "\n";
  std::cout << std::endl;
  std::cout << " Vincent.Richefeu@3sr-grenoble.fr AND Gael.Combe@3sr-grenoble.fr  " << "\n";
  std::cout << "\n\n";
  std::cout << " Version: " << TRACKER_VERSION << std::endl;

#if defined(_OPENMP)
  std::cout << " Accelerator: OPENMP" << std::endl;
#else
  std::cout << " Accelerator: None" << std::endl;
#endif

  std::cout << std::endl << std::endl;
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
      std::cout << std::endl;
      std::cout << "Help can be found in the folder 'examples' or 'doc'" << std::endl;
      std::cout << "  -v, --version   Display information about the soft" << std::endl;
      std::cout << "  -h, --help      This help" << std::endl;
      std::cout << "  -s, --sizes     Print size of common types (in bits) for this computer" << std::endl;
      std::cout << "  -g, --generate  Generate synthetic images for accuracy tests" << std::endl;
      std::cout << "  -d, --dialog    Use tracker to interactively extract information" << std::endl;
      std::cout << std::endl;
      exit(0);
    } else if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
      header();
      exit(0);
    } else if (strcmp(argv[1], "-s") == 0 || strcmp(argv[1], "--sizes") == 0) {
      std::cout << std::endl;
      std::cout << "unsigned char....." << sizeof(unsigned char) * 8 << " bits" << std::endl;
      std::cout << "unsigned short...." << sizeof(unsigned short) * 8 << " bits" << std::endl;
      std::cout << "unsigned int......" << sizeof(unsigned int) * 8 << " bits" << std::endl;
      std::cout << "float............." << sizeof(float) * 8 << " bits" << std::endl;
      std::cout << "double............" << sizeof(double) * 8 << " bits" << std::endl;
      std::cout << "uint16_t.........." << sizeof(uint16_t) * 8 << " bits" << std::endl;
      std::cout << "uint16_t..........." << sizeof(uint16_t) * 8 << " bits" << std::endl;
      std::cout << std::endl;
      exit(0);
    } else if (strcmp(argv[1], "-g") == 0 || strcmp(argv[1], "--generate") == 0) {
      int number_images = 0;
      double defx = 0.0, defy = 0.0, transx_pix = 0.0, transy_pix = 0.0, rot_deg = 0.0;
      double defx_inc = 0.0, defy_inc = 0.0, transx_pix_inc = 0.0, transy_pix_inc = 0.0, rot_deg_inc = 0.0;
      int w, h, octaves;
      double zoom, p;

      std::cout << std::endl;
      std::cout << "___Synthese of images." << std::endl;
      std::cout << "Number of images generated               : ";
      std::cin >> number_images;
      std::cout << "width (integer)                          : ";
      std::cin >> w;
      std::cout << "heigh (integer)                          : ";
      std::cin >> h;
      std::cout << std::endl;
      std::cout << "___Perlin noise. (octaves = 6, persistence = 0.5 and zoom = 5.0 seem to be ok)" << std::endl;
      std::cout << "octaves (integer)                        : ";
      std::cin >> octaves;
      std::cout << "persistence (double)                     : ";
      std::cin >> p;
      std::cout << "zoom (double)                            : ";
      std::cin >> zoom;
      std::cout << "___Transformations. (positive values correspond to left->right for x, top->bottom for y, "
                   "counterclockwise for rotation)"
                << std::endl;
      std::cout << "x-deformation increment (-, double) : ";
      std::cin >> defx_inc;
      std::cout << "y-deformation increment (-, double) : ";
      std::cin >> defy_inc;
      std::cout << "x-translation increment (pixels, double) : ";
      std::cin >> transx_pix_inc;
      std::cout << "y-translation increment (pixels, double) : ";
      std::cin >> transy_pix_inc;
      std::cout << "rotation increment (degrees, double)     : ";
      std::cin >> rot_deg_inc;

      // LOG in a file
      std::ofstream file("synth_log.txt");
      file << "___Synthese of images." << std::endl;
      file << "Number of images generated               : " << number_images << std::endl;
      file << "width (integer)                          : " << w << std::endl;
      file << "heigh (integer)                          : " << h << std::endl;
      file << std::endl;
      file << "___Perlin noise." << std::endl;
      file << "octaves (integer)                        : " << octaves << std::endl;
      file << "persistence (double)                     : " << p << std::endl;
      file << "zoom (double)                            : " << zoom << std::endl;
      file << "___Transformations." << std::endl;
      file << "x-deformation increment (pixels, double) : " << defx_inc << std::endl;
      file << "y-deformation increment (pixels, double) : " << defy_inc << std::endl;
      file << "x-translation increment (pixels, double) : " << transx_pix_inc << std::endl;
      file << "y-translation increment (pixels, double) : " << transy_pix_inc << std::endl;
      file << "rotation increment (degrees, double)     : " << rot_deg_inc << std::endl;

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
    TRACKER_LOG("Rotation are NOT TRACKED!" << std::endl);
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

  TRACKER_LOG("Number of threads: " << wanted_num_threads << std::endl);

  return 1;
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
  std::cout << std::endl;
  std::cout << " [ Enter Interactive Mode ]" << std::endl;
  std::cout << std::endl;
  std::cout << "Type 'help' for help" << std::endl;

  bool loaded = false;

  std::string token;
  std::cout << "tracker (interactive mode) > ";
  std::cin >> token;
  for (;;) {
    if (token == "help") {
      std::cout << "Type one of the following keywords:" << std::endl;
      std::cout << "  read_image  ..................... Read an image by using libtiff" << std::endl;
      std::cout << "  read_raw_image  ................. Read an image by using libraw" << std::endl;
      std::cout << "  DemosaicModel  .................. Set the value of DemosaicModel" << std::endl;
      std::cout << "  minmax  ......................... Show min and max values in the loaded image" << std::endl;
      std::cout << "  histo ........................... Compute the histogram of the loaded image" << std::endl;
      std::cout << "  image_check ..................... Create an image as it is 'seen' by tracker" << std::endl;
      // cout << "  draw_grains ..................... Superimpose grains over the thumb" << endl;
    } else if (token == "test") {
      im_index_current = 0;
      // read_raw_image(im_index_current, "../examples/ImageData/Disto_1.IIQ", true);
      read_image(im_index_current, "../examples/ImageData/TEST02.tiff", true);
      read_grains("../examples/ParticleImageTracking/dic_out_2.txt", true);
      image_div = 1;
      thumbnail tt(image, dimx, dimy, im_index_current, image_div);
      ColorTable ct;
      ct.setMinMax(0.0, 1.0);
      ct.savePpm("ColorScale.ppm");
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
      thumbnail tt(image, dimx, dimy, im_index_current, image_div);
      // tt.writePpm("visu.ppm", 0.3);
      tt.writeTiff("visuRaw.tiff", 0.0);
    } else if (token == "read_dic") {
      std::cout << "file path: ";
      std::string path;
      std::cin >> path;
      if (!grain.empty())
        grain.clear();
      read_grains(path.c_str(), true);
    } else if (token == "read_image") {
      std::cout << "file path: ";
      std::string path;
      std::cin >> path;
      im_index_current = 0;
      read_image(im_index_current, path.c_str(), true);
      if (dimx != 0)
        loaded = true;
    } else if (token == "read_raw_image") {
      std::cout << "file path: ";
      std::string path;
      std::cin >> path;
      im_index_current = 0;
      read_raw_image(im_index_current, path.c_str(), true);
      if (dimx != 0)
        loaded = true;
    } else if (token == "DemosaicModel") {
      std::cout << "value: ";
      std::cin >> DemosaicModel;
    } else if (token == "image_check") {
      if (!loaded) {
        std::cout << "Read an image first" << std::endl;
        token = "";
        continue;
      }
      im_index_current = 0;
      image_div = 1;
      thumbnail tt(image, dimx, dimy, im_index_current, image_div);
      tt.writeTiff("image_check.tiff", 0.0);
    } else if (token == "minmax") {
      if (!loaded) {
        std::cout << "Read an image first" << std::endl;
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
      std::cout << "minVal = " << minVal << std::endl;
      std::cout << "maxVal = " << maxVal << std::endl;
    } else if (token == "histo") {
      if (!loaded) {
        std::cout << "Read an image first" << std::endl;
        token = "";
        continue;
      }

      std::vector<size_t> hist(65535);
      for (size_t i = 0; i < 65535; i++) {
        hist[i] = 0;
      }

      for (int ix = 0; ix < dimx; ix++) {
        for (int iy = 0; iy < dimy; iy++) {
          hist[image[0][ix][iy]]++;
        }
      }

      std::ofstream file("histo.txt");
      for (size_t i = 0; i < 65535; i++) {
        file << hist[i] << std::endl;
      }
      std::cout << "Histogram saved in 'hist.txt'" << std::endl;
    } else if (token == "q" || token == "quit")
      break;

    std::cout << "tracker> ";
    std::cin >> token;
  }
  std::cout << " [ Exit Interactive Mode ]" << std::endl;
  exit(EXIT_SUCCESS);
}

/************************************************************************************************/
/*  				 PARTICLE TRACKING                                                                  */
/************************************************************************************************/

// The particle tracking procedure (default procedure)
void particle_tracking() {
  read_image(im_index_ref, iref, true); // Read reference image

  if (make_images) {
    create_image(iref);
  }

  double tbeg;
  int num_image;
  int igrain;
  int irescue;

  if (use_neighbour_list) {
    TRACKER_LOG("Build neighbour list ... " << std::flush);
    tbeg = get_time();
    for (int igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].neighbour.clear();
      find_neighbours(igrain);
    }
    TRACKER_LOG("[DONE in " << get_time() - tbeg << " seconds]" << std::endl);
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
        if (grain[igrain].masked) {
          continue;
        }
        follow_pattern_pixel(igrain);
#pragma omp critical
        {
          msg::loadbar(++progress, grain.size());
        }
      }

      std::cerr << std::endl;

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
        {
          msg::loadbar(++progress, num_to_be_rescued);
        }
      }
      std::cerr << std::endl;

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
      if (num_to_be_super_rescued > 0) {
        fprintf(stdout, "Try to super_rescue %d grains\n", num_to_be_super_rescued);
      }

      tbeg = get_time();

      progress = 0;
#pragma omp parallel for private(igrain) schedule(dynamic)
      for (irescue = 0; irescue < num_to_be_super_rescued; irescue++) {
        igrain = to_be_super_rescued[irescue];
        follow_pattern_super_rescue_pixel(igrain);
#pragma omp critical
        {
          msg::loadbar(++progress, num_to_be_rescued);
        }
      }
      std::cerr << std::endl;

      // Display
      for (irescue = 0; irescue < num_to_be_super_rescued; irescue++) {
        igrain = to_be_super_rescued[irescue];
        fprintf(stdout, "Grain %6d \t NCC = %6.4f  ->  ", igrain, grain[igrain].NCC);
        fprintf(stdout, "%6.4f \t diff = %6.4f \t ", grain[igrain].NCC_rescue,
                grain[igrain].NCC_rescue - grain[igrain].NCC);
        if (grain[igrain].NCC_rescue > NCC_min_super) {
          fprintf(stdout, "[\033[32mOK\033[0m]\n");
        } else {
          fprintf(stdout, "[\033[31mFAIL\033[0m]\n");
          fprintf(stdout, " Please, fix the pb manually\n");
        }
      }

      if (num_to_be_super_rescued > 0) {
        fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);
      }
    }

    // Sub-pixel
    if (subpixel) {
      fprintf(stdout, "Sub-pixel resolution ... \n");
      fflush(stdout);

      tbeg = get_time();

      progress = 0;
#pragma omp parallel for schedule(dynamic)
      for (igrain = 0; igrain < num_grains; igrain++) {
        if (grain[igrain].masked) {
          continue;
        }
        follow_pattern_subpixel_xyR(igrain);
#pragma omp critical
        {
          msg::loadbar(++progress, num_grains);
        }
      }
      std::cerr << std::endl;
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
      std::cout << std::endl;
      std::cout << "########" << std::endl << std::endl;
      std::cout << "Changing reference image from " << iref << " to " << num_image << std::endl << std::endl;
      std::cout << "########" << std::endl;

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
    std::cout << "@SolidMotionError, meanError = " << mean << ", stddevError = " << stddev << std::endl;
    std::cout << "angle  = " << angle << std::endl;
    std::cout << "xc     = " << xc << std::endl;
    std::cout << "yc     = " << yc << std::endl;
    std::cout << "xtrans = " << xtrans << std::endl;
    std::cout << "ytrans = " << ytrans << std::endl;
  }
}

// La nouvelle procédure discutée avec Gael pour que les corrélations rattée
// soient refaites mais cette fois avec une aide (une solution approchée) donnée
// par l'utiilsateur.
void particle_tracking_human_assisted() {
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
    } // Loop upix
  } // Loop rot

  grain[igrain].NCC = best_NCC;
  grain[igrain].NCC_rescue = best_NCC;

  if (best_NCC < NCC_min) {
#pragma omp critical
    {
      to_be_rescued[num_to_be_rescued++] = igrain;
    }
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
  std::vector<relative_coord_type> rescue_allowed(nbmax);

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
  } // Loop rot

  grain[igrain].NCC_rescue = best_NCC;

  if (best_NCC < NCC_min_super) {
#pragma omp critical
    {
      to_be_super_rescued[num_to_be_super_rescued++] = igrain;
    }
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
  std::vector<relative_coord_type> super_rescue_allowed(nbmax);

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

double NCC_to_minimize_xyR(std::vector<double> &X) {

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

  std::vector<double> interpol_values;

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
  std::vector<double> X(3); // Vector that hold the parameters to be optimized
  std::vector<double> DX(3);

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
  Powell<double(std::vector<double> &)> powell(NCC_to_minimize_xyR, subpix_tol);
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

int read_gmsh(const char *name) {
  std::ifstream meshFile(name);
  if (!meshFile) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  std::string token;
  char not_read[150];
  meshFile >> token;

  while (meshFile) {
    if (token == "$Nodes") {
      meshFile >> num_grains;
      TRACKER_SHOW(num_grains);

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
  std::ifstream grain_file(name);
  if (!grain_file) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  grain_file >> num_grains;
  TRACKER_SHOW(num_grains);

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

  std::cout << "Create grid with " << num_grains << " points" << std::endl;

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
  std::ofstream grain_file_out(name);
  if (!grain_file_out) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  if (simpleVersion) {
    grain_file_out << grain.size() << std::endl;
    for (size_t i = 0; i < grain.size(); i++) {
      grain_file_out << grain[i].refcoord_xpix + grain[i].dx << ' ' << grain[i].refcoord_ypix + grain[i].dy << ' '
                     << grain[i].refrot << ' ' << grain[i].radius_pix << '\n';
    }
  } else {
    grain_file_out << num_grains << std::endl;
    for (int i = 0; i < num_grains; i++) {
      grain_file_out << std::scientific << std::setprecision(15) << std::setw(6) << std::left
                     << grain[i].refcoord_xpix << ' '                               // 1
                     << std::setw(6) << std::left << grain[i].refcoord_ypix << ' '  // 2
                     << std::setw(23) << std::right << grain[i].refrot << ' '       // 3
                     << std::setw(23) << std::right << grain[i].radius_pix << ' '   // 4
                     << std::setw(23) << std::right << grain[i].dx << ' '           // 5
                     << std::setw(23) << std::right << grain[i].dy << ' '           // 6
                     << std::setw(23) << std::right << grain[i].drot << ' '         // 7
                     << std::setw(23) << std::right << grain[i].upix << ' '         // 8
                     << std::setw(23) << std::right << grain[i].vpix << ' '         // 9
                     << std::setw(23) << std::right << grain[i].rot_inc << ' '      // 10
                     << std::setw(23) << std::right << grain[i].NCC << ' '          // 11
                     << std::setw(23) << std::right << grain[i].NCC_rescue << ' '   // 12
                     << std::setw(23) << std::right << grain[i].NCC_subpix << '\n'; // 13
    }

    grain_file_out << std::endl;

    char fname[256];
    sprintf(fname, image_name, iref);
    grain_file_out << "# from image: " << fname << " --dateTime: " << imageData[im_index_ref].dateTime; //<< "\n";
    sprintf(fname, image_name, num);
    grain_file_out << "#   to image: " << fname << " --dateTime: " << imageData[im_index_current].dateTime; //<< "\n";
    // grain_file_out << fname << endl;
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
  std::ifstream pattern_file(name);
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
  std::ifstream command_file(name);
  if (!command_file) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  std::string token;
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
      std::string order;
      command_file >> order;
      if (order == "increasingRadius") {
        sort(grain.begin(), grain.end(), increasingRadius_order());
        save_grains("increasingRadius.txt", 0, true);
        std::cout << "Grains have been reordered with increasing radius in file increasingRadius.txt" << std::endl;
        exit(0);
      } else if (order == "decreasingRadius") {
        sort(grain.begin(), grain.end(), decreasingRadius_order());
        save_grains("decreasingRadius.txt", 0, true);
        std::cout << "Grains have been reordered with decreasing radius in file decreasingRadius.txt" << std::endl;
        exit(0);
      } else if (order == "increasingHeight") {
        sort(grain.begin(), grain.end(), increasingHeight_order());
        save_grains("increasingHeight.txt", 0, true);
        std::cout << "Grains have been reordered with increasing height in file increasingHeight.txt" << std::endl;
        exit(0);
      } else if (order == "mixed") {
        // std::random_shuffle(grain.begin(), grain.end());
        std::shuffle(grain.begin(), grain.end(), std::mt19937{std::random_device{}()});
        save_grains("mixed.txt", 0, true);
        std::cout << "Grains have been mixed in file mixed.txt" << std::endl;
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
      std::cout << "WARNING!! interp_mask_size is deprecated since version 0.5.0\n";
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
      std::string filename;
      command_file >> filename;
      std::string filename_to = "undist_" + filename;
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
      std::cout << "visu_mode = " << visu_mode << std::endl;
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
  std::ofstream command_file(name);
  if (!command_file) {
    fprintf(stderr, "Cannot open file %s\n", name);
    exit(0);
  }

  command_file << "image_name " << image_name << std::endl;
  command_file << "RawImages " << RawImages << std::endl;
  command_file << "DemosaicModel " << DemosaicModel << std::endl;

  command_file << grain_positions_command << std::endl; //// TODO  !!!!!!!!!!!!!!!!

  command_file << "iref " << iref << std::endl;
  command_file << "ibeg " << ibeg << std::endl;
  command_file << "iend " << iend << std::endl;
  command_file << "iinc " << iinc << std::endl;
  command_file << "iraz " << iraz << std::endl;

  command_file << "wanted_num_threads " << wanted_num_threads << std::endl;
  command_file << "procedure " << procedure << std::endl;
  command_file << "verbose_level " << verbose_level << std::endl;
  command_file << "rescue_level " << rescue_level << std::endl;
  command_file << "subpixel " << subpixel << std::endl;
  command_file << "rotations " << rotations << std::endl;

  command_file << "num_neighbour_max " << num_neighbour_max << std::endl;
  command_file << "neighbour_dist_pix " << neighbour_dist_pix << std::endl;
  command_file << "period_rebuild_neighbour_list " << period_rebuild_neighbour_list << std::endl;

  command_file << "image_interpolator "
               << "linear" << std::endl; ////// TODO !!!!!!!!!!!!!!
  command_file << "subpix_tol " << subpix_tol << std::endl;
  command_file << "initial_direction_dx " << initial_direction_dx << std::endl;
  command_file << "initial_direction_dy " << initial_direction_dy << std::endl;
  command_file << "initial_direction_dr " << initial_direction_dr << std::endl;

  command_file << "make_images " << make_images << std::endl;
  command_file << "image_size " << image_size << std::endl;
  command_file << "draw_angle " << draw_angle << std::endl;
  command_file << "draw_disp " << draw_disp << std::endl;
  command_file << "draw_rescued " << draw_rescued << std::endl;

  command_file << "search_zone.left " << search_zone.left << std::endl;
  command_file << "search_zone.right " << search_zone.right << std::endl;
  command_file << "search_zone.up " << search_zone.up << std::endl;
  command_file << "search_zone.down " << search_zone.down << std::endl;
  command_file << "search_zone.inc_rot " << search_zone.inc_rot << std::endl;
  command_file << "search_zone.num_rot " << search_zone.num_rot << std::endl;

  command_file << "NCC_min " << NCC_min << std::endl;
  command_file << "search_zone_rescue.left " << search_zone_rescue.left << std::endl;
  command_file << "search_zone_rescue.right " << search_zone_rescue.right << std::endl;
  command_file << "search_zone_rescue.up " << search_zone_rescue.up << std::endl;
  command_file << "search_zone_rescue.down " << search_zone_rescue.down << std::endl;
  command_file << "search_zone_rescue.inc_rot " << search_zone_rescue.inc_rot << std::endl;
  command_file << "search_zone_rescue.num_rot " << search_zone_rescue.num_rot << std::endl;

  command_file << "NCC_min_super " << NCC_min_super << std::endl;
  command_file << "search_zone_super_rescue.left " << search_zone_super_rescue.left << std::endl;
  command_file << "search_zone_super_rescue.right " << search_zone_super_rescue.right << std::endl;
  command_file << "search_zone_super_rescue.up " << search_zone_super_rescue.up << std::endl;
  command_file << "search_zone_super_rescue.down " << search_zone_super_rescue.down << std::endl;
  command_file << "search_zone_super_rescue.inc_rot " << search_zone_super_rescue.inc_rot << std::endl;
  command_file << "search_zone_super_rescue.num_rot " << search_zone_super_rescue.num_rot << std::endl;
  command_file << "use_neighbour_list " << use_neighbour_list << std::endl;

  command_file << pattern_command << std::endl; // TODO !!!!!!!!

  command_file << "subpix_center_threshold " << subpix_center_threshold << std::endl;
  command_file << "subpix_center_dr0 " << subpix_center_dr0 << std::endl;
  command_file << "subpix_center_rstep " << subpix_center_rstep << std::endl;

  // command_file << "image_numbers_corrDisto " << XXXX << endl;
  // command_file << "imposedDisplApprox_corrDisto " << XXXX << endl;

  command_file << "xc_corrDistor " << disto_parameters[0] << std::endl;
  command_file << "yc_corrDistor " << disto_parameters[1] << std::endl;
  command_file << "K1_corrDistor " << disto_parameters[2] << std::endl;
  command_file << "K2_corrDistor " << disto_parameters[3] << std::endl;
  command_file << "K3_corrDistor " << disto_parameters[4] << std::endl;
  command_file << "P1_corrDistor " << disto_parameters[5] << std::endl;
  command_file << "P2_corrDistor " << disto_parameters[6] << std::endl;
  command_file << "P3_corrDistor " << disto_parameters[7] << std::endl;

  command_file << "dxc_corrDistor " << disto_parameters_perturb[0] << std::endl;
  command_file << "dyc_corrDistor " << disto_parameters_perturb[1] << std::endl;
  command_file << "dK1_corrDistor " << disto_parameters_perturb[2] << std::endl;
  command_file << "dK2_corrDistor " << disto_parameters_perturb[3] << std::endl;
  command_file << "dK3_corrDistor " << disto_parameters_perturb[4] << std::endl;
  command_file << "dP1_corrDistor " << disto_parameters_perturb[5] << std::endl;
  command_file << "dP2_corrDistor " << disto_parameters_perturb[6] << std::endl;
  command_file << "dP3_corrDistor " << disto_parameters_perturb[7] << std::endl;

  command_file << "! END OF FILE";
}

// #include "image_io.cpp"

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
    std::cerr << "Cannot read tiff file named '" << name << "'" << std::endl;
    exit(0);
  }
  uint32_t w, h;
  size_t npixels;
  uint32_t *raster;

  char *infobuf;
  if (TIFFGetField(tif, TIFFTAG_DATETIME, &infobuf))
    imageData[i].dateTime = std::string(infobuf);
  else
    imageData[i].dateTime = "dateTime unknown";

  imageData[i].iso_speed = 0.0;
  imageData[i].shutter = 0.0;
  imageData[i].aperture = 0.0;
  imageData[i].focal_len = 0.0;
  imageData[i].shot_order = 0;

  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  std::cout << w << "x" << h << std::endl;
  npixels = w * h;

  raster = (uint32_t *)_TIFFmalloc(npixels * sizeof(uint32_t));
  if (raster != NULL) {
    if (!TIFFReadRGBAImage(tif, w, h, raster, 0)) {
      std::cerr << "Cannot read tiff file named '" << name << "'" << std::endl;
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
      if (r < 0.0) {
        r = 0.0;
        std::cerr << "@read_image, r < 0.0\n";
      }
      if (r > 1.0) {
        r = 1.0;
        std::cerr << "@read_image, r > 1.0\n";
      }
      if (g < 0.0) {
        g = 0.0;
        std::cerr << "@read_image, g < 0.0\n";
      }
      if (g > 1.0) {
        g = 1.0;
        std::cerr << "@read_image, g > 1.0\n";
      }
      if (b < 0.0) {
        b = 0.0;
        std::cerr << "@read_image, b < 0.0\n";
      }
      if (b > 1.0) {
        b = 1.0;
        std::cerr << "@read_image, b > 1.0\n";
      }
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


// il faut voir ici http://www.libraw.org/node/555 pour voir comment proceder
// ici aussi :
// http://stackoverflow.com/questions/22355491/libraw-is-making-my-images-too-bright-compared-to-nikons-own-converter
void read_raw_image(int i, const char *name, bool first_time) {

  LibRaw iProcessor;
  int IO_error = iProcessor.open_file(name);
  if (IO_error != 0) {
    std::cout << std::endl;
    std::cout << "The image named " << name << " does not exist..." << std::endl;
    std::cout << "   ## PROGRAM STOPPED ##" << std::endl;
    exit(1);
  }

  double tbeg = get_time();
  fprintf(stdout, "Read and demosaicing image named %s (DemosaicModel = %d)... \n", name, DemosaicModel);

  if (first_time) {
    dimx = iProcessor.imgdata.sizes.width;
    dimy = iProcessor.imgdata.sizes.height;
    std::cout << "Size = " << dimx << "x" << dimy << std::endl;

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

  auto timestamp2string = [](time_t rawtime) -> std::string {
    struct tm timeinfo;
    localtime_r(&rawtime, &timeinfo); // Version thread-safe (POSIX)
    std::stringstream ss;
    ss << std::put_time(&timeinfo, "%a %b %d %H:%M:%S %Y\n");
    return ss.str();
  };
  imageData[i].dateTime = timestamp2string(iProcessor.imgdata.other.timestamp);

  iProcessor.unpack();
  uint16_t MinGray = 65535, MaxGray = 0;
  double fact = 1.0 / (double)(iProcessor.imgdata.color.maximum);
  // std::cout << "iProcessor.imgdata.color.maximum = " << iProcessor.imgdata.color.maximum << '\n';
  // std::cout << "fact = " << fact << '\n';

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

}

// This function has to be rewritten for a pgm (or tiff) output
void undistor_image(const char * /*name_from*/, const char * /*name_to*/) {
  read_image(im_index_ref, grid_image_name.c_str(), true);

  double xu, yu;
  for (int y = 0; y < dimy; ++y) {
    for (int x = 0; x < dimx; ++x) {
      undistor(&disto_parameters[0], x, y, xu, yu);
      xu = nearest(xu);
      yu = nearest(yu);
      if (xu < 0. || xu >= dimx || yu < 0. || yu >= dimy)
        continue;
      image[im_index_current][(int)xu][(int)yu] = image[im_index_ref][x][y];
    }
  }

  create_image_Netbpm(0);
}



void draw_circle(std::vector<std::vector<uint16_t>> &imThumb, int x_centre, int y_centre, int radius, uint16_t col) {
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
  std::cout << "Thumb size: " << imW << "x" << imH << std::endl;

  std::vector<std::vector<uint16_t>> imThumb;
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
  std::cout << "min_grey = " << min_grey << std::endl;
  std::cout << "max_grey = " << max_grey << std::endl;
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
  std::cout << "Create image " << num << " for check ... " << std::flush;
  char name[256];
  sprintf(name, "check_%d.pgm", num);
  std::ofstream file(name, std::ios::binary);
  file << "P5\n";
  file << imW << " " << imH << std::endl;
  file << "255" << std::endl;

  double conv = 255.0 / 65535.0;

  for (int y = 0; y < imH; y++) {
    for (int x = 0; x < imW; x++) {
      unsigned char v = nearest(conv * (imThumb[x][y]));
      file.write((const char *)&v, sizeof(unsigned char));
    }
  }

  std::cout << "done." << std::endl;
}

void create_image(int num) { create_image_Netbpm(num); }

/************************************************************************************************/
/*                                CORRECTION OF DISTORSION                                      */
/************************************************************************************************/

// http://en.wikipedia.org/wiki/Distortion_(optics)
// + pdf file named 'distortion' in the folder 'refs'
void distor(double *X, const double xu, const double yu, double &xd, double &yd) {
  double Xc = X[0];
  double Yc = X[1];
  double K1 = X[2];
  double K2 = X[3];
  double K3 = X[4];
  double P1 = X[5];
  double P2 = X[6];
  double P3 = X[7];

  double dx = xu - Xc;
  double dy = yu - Yc;
  double r2 = dx * dx + dy * dy;
  double r4 = r2 * r2;
  double r6 = r4 * r2;

  // Brown's lens distortion formula
  // http://www.ise.pw.edu.pl/~rrom/SPIE/SPIE8903-WILGA2013/8903/CD/WI13/DATA/8903_92.PDF
  xd = Xc + dx * (1.0 + K1 * r2 + K2 * r4 + K3 * r6) +
       (P1 * (r2 + 2.0 * dx * dx) + 2.0 * P2 * dx * dy) * (1.0 + P3 * r2);
  yd = Yc + dy * (1.0 + K1 * r2 + K2 * r4 + K3 * r6) +
       (P1 * (r2 + 2.0 * dy * dy) + 2.0 * P2 * dx * dy) * (1.0 + P3 * r2);
}

void diff_distor(double *X, const double xu, const double yu, double &dxd_dxu, double &dyd_dyu) {
  double Xc = X[0];
  double Yc = X[1];
  double K1 = X[2];
  double K2 = X[3];
  double K3 = X[4];
  double P1 = X[5];
  double P2 = X[6];
  double P3 = X[7];

  double dx = xu - Xc;
  double dy = yu - Yc;
  double r2 = dx * dx + dy * dy;
  double r4 = r2 * r2;
  double r6 = r4 * r2;

  dxd_dxu = (1.0 + K1 * r2 + K2 * r4 + K3 * r6) + 2.0 * dx * dx * (2.0 * K2 * r2 + K1 + 3 * K3 * r4) +
            (6.0 * P1 * dx + 2.0 * P2 * dy) * (1.0 + P3 * r2) +
            (P1 * (r2 + 2.0 * dx * dx) + 2.0 * P2 * dx * dy) * P3 * 2.0 * dx;

  dyd_dyu = (1.0 + K1 * r2 + K2 * r4 + K3 * r6) + 2.0 * dy * dy * (2.0 * K2 * r2 + K1 + 3 * K3 * r4) +
            (6.0 * P1 * dy + 2.0 * P2 * dx) * (1.0 + P3 * r2) +
            (P1 * (r2 + 2.0 * dy * dy) + 2.0 * P2 * dy * dx) * P3 * 2.0 * dy;
}

int undistor(double *X, const double xd, const double yd, double &xu, double &yu, int niterMax, double tol,
             int nconv) {
  if (fake_undistor == 1) {
    // Fonction de distorsion utilisée à l'envers
    // It is how undistor is used most often. And it is okay! not that sure-gael
    distor(&X[0], xd, yd, xu, yu);
    return 1;
  }

  // initial guess
  xu = xd;
  yu = yd;

  // convergence criterion
  int kconv = 0;
  int nbiter = 0;

  double fx, fy;
  double dxd_dxu, dyd_dyu;
  double deltax, deltay;

  for (int i = 0; i < niterMax; i++) {
    nbiter++;
    // derivatives
    diff_distor(&X[0], xu, yu, dxd_dxu, dyd_dyu);
    // re-estimate
    distor(&X[0], xu, yu, fx, fy);
    if (fabs(dxd_dxu) > 1e-15)
      deltax = (xd - fx) / dxd_dxu;
    else
      deltax = 0.0;
    if (fabs(dyd_dyu) > 1e-15)
      deltay = (yd - fy) / dyd_dyu;
    else
      deltay = 0.0;
    xu += deltax;
    yu += deltay;

    // check convergence
    if (deltax < tol && deltay < tol)
      kconv++;
    else
      kconv = 0;
    if (kconv >= nconv)
      break;
  }
  return nbiter;
}

// Cette fonction construit des pairs de points sur une grille
// qui sont distant de plus de equiProj_dist_min et de moins de equiProj_dist_max
void precompute_paires() {
  List_Grains_i.clear();
  List_Grains_j.clear();

  for (int igrain = 0; igrain < num_grains; igrain++) {
    double xi = grain[igrain].refcoord_xpix;
    double yi = grain[igrain].refcoord_ypix;
    for (int jgrain = igrain + 1; jgrain < num_grains; jgrain++) {
      double xj = grain[jgrain].refcoord_xpix;
      double yj = grain[jgrain].refcoord_ypix;
      double xij = (xj - xi);
      double yij = (yj - yi);
      double dij = sqrt(xij * xij + yij * yij);
      if (dij >= equiProj_dist_min && dij <= equiProj_dist_max) {
        List_Grains_i.push_back(igrain);
        List_Grains_j.push_back(jgrain);
      }
    }
  }

  std::cout << "Number of points for equiprojectivity computations: " << List_Grains_i.size() << '\n';
}

// La fonction a minimiser si on veut utiliser le critère d'équiprojectivité.
// Pour il faut que les déplacements entre les image correspondent réellement à un mouvement
// de solide rigide (déplacement d'une grande plaque).
//
// Remarque suite aux tests de Gael :
// La minimisation de sum(|err|) se comporte mal! (solution très sensible aux perturbations de départ)
// Il vaut mieux minimiser sum(err). Mais on préfère finalement minimiser max(|err|)
double disto_to_minimize_Equiproj(std::vector<double> &X) {
  int igrain, jgrain;
  double A1x, A1y, B1x, B1y;
  double A2x, A2y, B2x, B2y;
  double AAx, AAy, BBx, BBy, ABx, ABy;
  double func = 0.0;
  double val = 0.0;

  int nbval = 0;
  for (size_t i_image = 1; i_image < image_numbers_corrDisto.size(); i_image++) {
    for (size_t il = 0; il < List_Grains_i.size(); il++) {
      igrain = List_Grains_i[il];
      jgrain = List_Grains_j[il];
      if (NCC_subpix_corrDisto[i_image][igrain] < NCC_min || NCC_subpix_corrDisto[i_image][jgrain] < NCC_min)
        continue;

      nbval++;

      // undistor(&X[0], (double)(grain[igrain].refcoord_xpix), (double)(grain[igrain].refcoord_ypix), A1x, A1y);
      A1x = (double)(grain[igrain].refcoord_xpix);
      A1y = (double)(grain[igrain].refcoord_ypix);
      undistor(&X[0], (double)(grain[igrain].refcoord_xpix + dx_corrDisto[i_image][igrain]),
               (double)(grain[igrain].refcoord_ypix + dy_corrDisto[i_image][igrain]), A2x, A2y);
      // undistor(&X[0], (double)(grain[jgrain].refcoord_xpix), (double)(grain[jgrain].refcoord_ypix), B1x, B1y);
      B1x = (double)(grain[jgrain].refcoord_xpix);
      B1y = (double)(grain[jgrain].refcoord_ypix);
      undistor(&X[0], (double)(grain[jgrain].refcoord_xpix + dx_corrDisto[i_image][jgrain]),
               (double)(grain[jgrain].refcoord_ypix + dy_corrDisto[i_image][jgrain]), B2x, B2y);

      ABx = B1x - A1x;
      ABy = B1y - A1y;
      double AB = sqrt(ABx * ABx + ABy * ABy);
      ABx /= AB;
      ABy /= AB;

      BBx = B2x - B1x;
      BBy = B2y - B1y;
      AAx = A2x - A1x;
      AAy = A2y - A1y;

      // version sum(err)

      // on prend seulement la plus grande erreur en valeur absolue
      val = (ABx * (BBx - AAx)) + (ABy * (BBy - AAy));
      func += val * val;
      // if (val > func) func = val;
    }
  }
  // return fabs(func);
  return (func);
}

double Equiproj_Error_Computation(std::vector<double> &X) {
  int igrain, jgrain;
  double A1x, A1y, B1x, B1y;
  double A2x, A2y, B2x, B2y;
  double AAx, AAy, BBx, BBy, ABx, ABy;
  double func = 0.0;

  for (size_t i_image = 1; i_image < image_numbers_corrDisto.size(); i_image++) {
    for (size_t il = 0; il < List_Grains_i.size(); il++) {
      igrain = List_Grains_i[il];
      jgrain = List_Grains_j[il];
      if (NCC_subpix_corrDisto[i_image][igrain] < NCC_min || NCC_subpix_corrDisto[i_image][jgrain] < NCC_min)
        continue;

      undistor(&X[0], (double)(grain[igrain].refcoord_xpix), (double)(grain[igrain].refcoord_ypix), A1x, A1y);
      undistor(&X[0], (double)(grain[igrain].refcoord_xpix + dx_corrDisto[i_image][igrain]),
               (double)(grain[igrain].refcoord_ypix + dy_corrDisto[i_image][igrain]), A2x, A2y);
      undistor(&X[0], (double)(grain[jgrain].refcoord_xpix), (double)(grain[jgrain].refcoord_ypix), B1x, B1y);
      undistor(&X[0], (double)(grain[jgrain].refcoord_xpix + dx_corrDisto[i_image][jgrain]),
               (double)(grain[jgrain].refcoord_ypix + dy_corrDisto[i_image][jgrain]), B2x, B2y);

      ABx = B1x - A1x;
      ABy = B1y - A1y;
      BBx = B2x - B1x;
      BBy = B2y - B1y;
      AAx = A2x - A1x;
      AAy = A2y - A1y;

      // version sum(err)
      func += abs((ABx * (BBx - AAx)) + (ABy * (BBy - AAy)));
    }
  }
  return func;
}

void correction_distortion() {
  std::cout << "CORRECTION OF DISTORTION" << std::endl;

  // Reserve memory
  dx_corrDisto.resize(image_numbers_corrDisto.size());
  dy_corrDisto.resize(image_numbers_corrDisto.size());
  NCC_subpix_corrDisto.resize(image_numbers_corrDisto.size());
  for (size_t i = 0; i < image_numbers_corrDisto.size(); i++) {
    dx_corrDisto[i].resize(num_grains);
    dy_corrDisto[i].resize(num_grains);
    NCC_subpix_corrDisto[i].resize(num_grains);
  }

  int igrain;
  double index_ini, index_end, index_gain;

  std::ofstream logfile("corrDisto.log");

  read_image(0, image_numbers_corrDisto[0], true); // Read reference image
  iref = image_numbers_corrDisto[0];
  do_precomputations();
  precompute_paires();

  if (make_images)
    create_image(0);

  std::cout << "Tracking displacements... " << std::endl << std::flush;
  double tbeg;
  for (size_t i_image = 1; i_image < image_numbers_corrDisto.size(); i_image++) {

    read_image(1, image_numbers_corrDisto[i_image]);

    // Accounting for the estimated displacement
    for (igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].dx = imposed_displ_x_pix[i_image];
      grain[igrain].dy = imposed_displ_y_pix[i_image];
    }

    // ======== Performing a DIC
    fprintf(stdout, "Follow %d grains ... (from image number %d to number %d)\n", num_grains, iref,
            image_numbers_corrDisto[i_image]);
    fflush(stdout);
    tbeg = get_time();

    progress = 0;
#pragma omp parallel for
    for (igrain = 0; igrain < num_grains; igrain++) {
      grain[igrain].reset();
      follow_pattern_pixel(igrain);
#pragma omp critical
      {
        msg::loadbar(++progress, grain.size());
      }
    }
    std::cerr << std::endl;
    fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);

    // NO RESCUE IS ATTEMPTED !!?

    // ======== Sub-pixel
    fprintf(stdout, "Sub-pixel resolution ... \n");
    fflush(stdout);
    tbeg = get_time();

    progress = 0;
#pragma omp parallel for
    for (igrain = 0; igrain < num_grains; igrain++) {
      follow_pattern_subpixel_xyR(igrain);
#pragma omp critical
      {
        msg::loadbar(++progress, grain.size());
      }
    }
    std::cerr << std::endl;
    fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);

    int nbFailed = 0;
    for (igrain = 0; igrain < num_grains; igrain++) {
      dx_corrDisto[i_image][igrain] = grain[igrain].dx;
      dy_corrDisto[i_image][igrain] = grain[igrain].dy;
      NCC_subpix_corrDisto[i_image][igrain] = grain[igrain].NCC_subpix;
      if (grain[igrain].NCC_subpix < NCC_min)
        nbFailed++;
    }
    std::cout << "Failure rate = " << 100.0 * (double)nbFailed / (double)num_grains << "%" << std::endl;

    // save_grains(i_image);
    save_grains(image_numbers_corrDisto[i_image]);

  } // for loop i_image
  std::cout << "done." << std::endl;

  // The parameters (initial guess) and perturbations are set by the user
  if (disto_parameters[0] < 0.0)
    disto_parameters[0] = 0.5 * (double)dimx;
  if (disto_parameters[1] < 0.0)
    disto_parameters[1] = 0.5 * (double)dimy;

  index_ini = Equiproj_Error_Computation(disto_parameters);

  // Minimization
  int npairs = List_Grains_i.size();
  std::cout << "Numbers of pairs used to undistor : " << npairs << std::endl;

  std::cout << "Minimizing equiprojectivity criterion... " << std::flush;
  Powell<double(std::vector<double> &)> powell(disto_to_minimize_Equiproj, 1e-12);
  disto_parameters = powell.minimize(disto_parameters, disto_parameters_perturb);
  std::cout << "done." << std::endl << std::endl;

  index_end = Equiproj_Error_Computation(disto_parameters);
  index_gain = index_ini / index_end - 1.0;
  logfile << std::scientific << std::setprecision(10);
  logfile << "fake_undistor            " << fake_undistor << std::endl;
  logfile << "Initial distorsion index " << index_ini << std::endl;
  logfile << "Final distorsion index   " << index_end << std::endl;
  logfile << "Index gain (>0 is good)  " << index_gain << std::endl;
  logfile << "xc_corrDistor            " << disto_parameters[0] << std::endl;
  logfile << "yc_corrDistor            " << disto_parameters[1] << std::endl;
  logfile << "K1_corrDistor            " << disto_parameters[2] << std::endl;
  logfile << "K2_corrDistor            " << disto_parameters[3] << std::endl;
  logfile << "K3_corrDistor            " << disto_parameters[4] << std::endl;
  logfile << "P1_corrDistor            " << disto_parameters[5] << std::endl;
  logfile << "P2_corrDistor            " << disto_parameters[6] << std::endl;
  logfile << "P3_corrDistor            " << disto_parameters[7] << std::endl;

  std::cout << std::scientific << std::setprecision(10);
  std::cout << "fake_undistor            " << fake_undistor << std::endl;
  std::cout << "Initial distorsion index " << index_ini << std::endl;
  std::cout << "Final distorsion index   " << index_end << std::endl;
  std::cout << "Index gain (>0 is good)  " << index_gain << std::endl;
  std::cout << "xc_corrDistor            " << disto_parameters[0] << std::endl;
  std::cout << "yc_corrDistor            " << disto_parameters[1] << std::endl;
  std::cout << "K1_corrDistor            " << disto_parameters[2] << std::endl;
  std::cout << "K2_corrDistor            " << disto_parameters[3] << std::endl;
  std::cout << "K3_corrDistor            " << disto_parameters[4] << std::endl;
  std::cout << "P1_corrDistor            " << disto_parameters[5] << std::endl;
  std::cout << "P2_corrDistor            " << disto_parameters[6] << std::endl;
  std::cout << "P3_corrDistor            " << disto_parameters[7] << std::endl;

  std::ofstream errorFile("errors.dat");
  double xerrormax = -1e20;
  double yerrormax = -1e20;
  double xd, yd, xu, yu;
  double dxd = (double)dimx / 1000.0;
  double dyd = (double)dimy / 1000.0;
  if (dxd < 1.0)
    dxd = 1.0;
  if (dyd < 1.0)
    dyd = 1.0;
  for (xd = 0; xd < dimx; xd += dxd) {
    for (yd = 0; yd < dimy; yd += dyd) {
      undistor(&disto_parameters[0], xd, yd, xu, yu);
      double xerr = xu - xd;
      double yerr = yu - yd;
      if (fabs(xerr) > xerrormax)
        xerrormax = fabs(xerr);
      if (fabs(yerr) > yerrormax)
        yerrormax = fabs(yerr);
      errorFile << xd << " " << yd << " " << xerr << " " << yerr << std::endl;
    }
    errorFile << std::endl;
  }
  logfile << "xerrormax                " << xerrormax << std::endl;
  logfile << "yerrormax                " << yerrormax << std::endl;
  std::cout << "xerrormax                " << xerrormax << std::endl;
  std::cout << "yerrormax                " << yerrormax << std::endl;

  std::ofstream errorBoxFile("errorbox.dat");
  yd = 0;
  for (xd = 0; xd < dimx; xd += dxd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  xd = dimx - 1;
  for (yd = 0; yd < dimy; yd += dyd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  yd = dimy - 1;
  for (xd = dimx - 1; xd >= 0; xd -= dxd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  xd = 0;
  for (yd = dimy - 1; yd >= 0; yd -= dyd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }

  xd = yd = 0;
  undistor(&disto_parameters[0], xd, yd, xu, yu);
  errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
}

void correction_distortion_grid() {
  std::cout << "CORRECTION OF DISTORSION -- GRID METHOD" << std::endl;
  read_image(im_index_ref, grid_image_name.c_str(), true); // Read grid image (to get dimx and dimy)
  std::ofstream logfile("corrDisto.log");
  double index = disto_to_minimize_grid(disto_parameters);
  logfile << "Initial distorsion index " << index << std::endl;
  std::cout << "Initial distorsion index " << index << std::endl;

  // The parameters (initial guess) and perturbations are set by the user
  if (disto_parameters[0] < 0.0)
    disto_parameters[0] = 0.5 * (double)dimx;
  if (disto_parameters[1] < 0.0)
    disto_parameters[1] = 0.5 * (double)dimy;

  // Minimization
  Powell<double(std::vector<double> &)> powell(disto_to_minimize_grid, 1e-8);
  disto_parameters = powell.minimize(disto_parameters, disto_parameters_perturb);

  index = disto_to_minimize_grid(disto_parameters);
  logfile << "Final distorsion index   " << index << std::endl;
  logfile << "xc_corrDistor            " << disto_parameters[0] << std::endl;
  logfile << "yc_corrDistor            " << disto_parameters[1] << std::endl;
  logfile << "K1_corrDistor            " << disto_parameters[2] << std::endl;
  logfile << "K2_corrDistor            " << disto_parameters[3] << std::endl;
  logfile << "K3_corrDistor            " << disto_parameters[4] << std::endl;
  logfile << "P1_corrDistor            " << disto_parameters[5] << std::endl;
  logfile << "P2_corrDistor            " << disto_parameters[6] << std::endl;
  logfile << "P3_corrDistor            " << disto_parameters[7] << std::endl;

  std::cout << "Final distorsion index   " << index << std::endl;
  std::cout << "xc_corrDistor            " << disto_parameters[0] << std::endl;
  std::cout << "yc_corrDistor            " << disto_parameters[1] << std::endl;
  std::cout << "K1_corrDistor            " << disto_parameters[2] << std::endl;
  std::cout << "K2_corrDistor            " << disto_parameters[3] << std::endl;
  std::cout << "K3_corrDistor            " << disto_parameters[4] << std::endl;
  std::cout << "P1_corrDistor            " << disto_parameters[5] << std::endl;
  std::cout << "P2_corrDistor            " << disto_parameters[6] << std::endl;
  std::cout << "P3_corrDistor            " << disto_parameters[7] << std::endl;

  std::ofstream gridFile("grid.dat");
  double xd, yd, xu, yu;
  for (int iy = 0; iy < ny_grid_disto; iy++) {
    for (int ix = 0; ix < nx_grid_disto; ix++) {
      int igrain = iy * nx_grid_disto + ix;
      xd = grain[igrain].refcoord_xpix;
      yd = grain[igrain].refcoord_ypix;
      undistor(&disto_parameters[0], xd, yd, xu, yu);
      gridFile << xd << " " << yd << " " << xu << " " << yu << std::endl;
    }
    gridFile << std::endl;
  }
  gridFile << std::endl;
  for (int ix = 0; ix < nx_grid_disto; ix++) {
    for (int iy = 0; iy < ny_grid_disto; iy++) {
      int igrain = iy * nx_grid_disto + ix;
      xd = grain[igrain].refcoord_xpix;
      yd = grain[igrain].refcoord_ypix;
      undistor(&disto_parameters[0], xd, yd, xu, yu);
      gridFile << xd << " " << yd << " " << xu << " " << yu << std::endl;
    }
    gridFile << std::endl;
  }

  std::ofstream errorFile("errors.dat");
  double dxd = (double)dimx / 100.0;
  double dyd = (double)dimy / 100.0;
  if (dxd < 1.0)
    dxd = 1.0;
  if (dyd < 1.0)
    dyd = 1.0;
  for (xd = 0; xd < dimx; xd += dxd) {
    for (yd = 0; yd < dimy; yd += dyd) {
      undistor(&disto_parameters[0], xd, yd, xu, yu);
      errorFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
    }
    errorFile << std::endl;
  }

  std::ofstream errorBoxFile("errorbox.dat");
  yd = 0;
  for (xd = 0; xd < dimx; xd += dxd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  xd = dimx - 1;
  for (yd = 0; yd < dimy; yd += dyd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  yd = dimy - 1;
  for (xd = dimx - 1; xd >= 0; xd -= dxd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }
  xd = 0;
  for (yd = dimy - 1; yd >= 0; yd -= dyd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;
  }

  xd = yd = 0;
  undistor(&disto_parameters[0], xd, yd, xu, yu);
  errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << std::endl;

  // show_distor();
  undistor_image(grid_image_name.c_str(), "undisto.png");
}

// To force horizontal lines to be horizontal, and vertical lines to be vertical
double disto_to_minimize_grid(std::vector<double> &X) {
  double sum = 0.0;
  int n = 0;
  double x0, y0, xu, yu;
  int igrain;
  for (int iy = 0; iy < ny_grid_disto; iy++) {
    igrain = iy * nx_grid_disto;
    undistor(&X[0], (double)(grain[igrain].refcoord_xpix), (double)(grain[igrain].refcoord_ypix), x0, y0);
    for (int ix = 0; ix < nx_grid_disto; ix++) {
      igrain = iy * nx_grid_disto + ix;
      undistor(&X[0], (double)(grain[igrain].refcoord_xpix), (double)(grain[igrain].refcoord_ypix), xu, yu);
      sum += (yu - y0) * (yu - y0);
      n++;
    }
  }
  for (int ix = 0; ix < nx_grid_disto; ix++) {
    undistor(&X[0], (double)(grain[ix].refcoord_xpix), (double)(grain[ix].refcoord_ypix), x0, y0);
    for (int iy = 0; iy < ny_grid_disto; iy++) {
      igrain = iy * nx_grid_disto + ix;
      undistor(&X[0], (double)(grain[igrain].refcoord_xpix), (double)(grain[igrain].refcoord_ypix), xu, yu);
      sum += (xu - x0) * (xu - x0);
      n++;
    }
  }
  return sum / (double)n;
}

void setSolidSolutionFromRef(double DX, double DY, double DROT) {
  for (int igrain = 0; igrain < num_grains; igrain++) {
    grain[igrain].dx = DX;
    grain[igrain].dy = DY;
    grain[igrain].drot = DROT;
  }
}

double SolidMotionError_to_minimize(std::vector<double> &X) {
  double angle = X[0];
  double xc = X[1];
  double yc = X[2];
  double xtrans = X[3];
  double ytrans = X[4];

  double c = cos(angle);
  double s = sin(angle);

  double Rxx = c;
  double Rxy = s;
  double Ryx = -s;
  double Ryy = c;

  double error = 0.0;
  for (int igrain = 0; igrain < num_grains; igrain++) {
    double deltax = (double)grain[igrain].refcoord_xpix - xc;
    double deltay = (double)grain[igrain].refcoord_ypix - yc;
    double xt = Rxx * deltax + Rxy * deltay + xtrans;
    double yt = Ryx * deltax + Ryy * deltay + ytrans;

    double dx = xt - (grain[igrain].refcoord_xpix + grain[igrain].dx);
    double dy = yt - (grain[igrain].refcoord_ypix + grain[igrain].dy);
    error += dx * dx + dy * dy;
  }
  return error /= (double)num_grains;
}

void SolidMotionError(double &angle, double &xc, double &yc, double &xtrans, double &ytrans, double &mean,
                      double &stddev) {
  std::vector<double> X(5);
  X[0] = angle;
  X[1] = xc;
  X[2] = yc;
  X[3] = xtrans;
  X[4] = ytrans;

  std::vector<double> dX(5);
  dX[0] = 0.0; // 1.0e-6;
  dX[1] = 0.0; // 1.0e-6;
  dX[2] = 0.0; // 1.0e-6;
  dX[3] = 1.0e-6;
  dX[4] = 1.0e-6;

  // Minimization
  Powell<double(std::vector<double> &)> powell(SolidMotionError_to_minimize, 1e-8);
  X = powell.minimize(X, dX);

  angle = X[0];
  xc = X[1];
  yc = X[2];
  xtrans = X[3];
  ytrans = X[4];

  // Calcul des indicateurs
  double c = cos(angle);
  double s = sin(angle);

  double Rxx = c;
  double Rxy = s;
  double Ryx = -s;
  double Ryy = c;

  mean = 0.0;
  for (int igrain = 0; igrain < num_grains; igrain++) {
    double deltax = (double)grain[igrain].refcoord_xpix - xc;
    double deltay = (double)grain[igrain].refcoord_ypix - yc;
    double xt = Rxx * deltax + Rxy * deltay + xtrans;
    double yt = Ryx * deltax + Ryy * deltay + ytrans;

    double dx = xt - (grain[igrain].refcoord_xpix + grain[igrain].dx);
    double dy = yt - (grain[igrain].refcoord_ypix + grain[igrain].dy);
    mean += dx * dx + dy * dy;
  }
  mean /= (double)num_grains;

  stddev = 0.0;
  for (int igrain = 0; igrain < num_grains; igrain++) {
    double deltax = (double)grain[igrain].refcoord_xpix - X[1];
    double deltay = (double)grain[igrain].refcoord_ypix - X[2];
    double xt = Rxx * deltax + Rxy * deltay + X[3];
    double yt = Ryx * deltax + Ryy * deltay + X[4];

    double dx = xt - (grain[igrain].refcoord_xpix + grain[igrain].dx);
    double dy = yt - (grain[igrain].refcoord_ypix + grain[igrain].dy);
    stddev += (dx * dx + dy * dy) - mean;
  }
  stddev /= (double)num_grains;
}

/************************************************************************************************/
/*                                     GRAY LEVEL ANALYSIS                                      */
/************************************************************************************************/

//#include "moment.hpp"



void moment(std::vector<double> &data, data_stat &stat) {
  size_t n = data.size();

  if (n <= 1)
    std::cerr << "@moment, number of data too small (min = 2)\n";
  double s = 0.0;
  for (size_t j = 0; j < n; j++)
    s += data[j];
  stat.ave = s / (double)n;
  stat.adev = stat.var = stat.skew = stat.curt = 0.0;
  double ep = 0.0;
  double p = 0.0;
  for (size_t j = 0; j < n; j++) {
    stat.adev += fabs(s = data[j] - stat.ave);
    ep += s;
    stat.var += (p = s * s);
    stat.skew += (p *= s);
    stat.curt += (p *= s);
  }
  stat.adev /= (double)n;
  stat.var = (stat.var - ep * ep / (double)n) / (double)(n - 1);
  stat.sdev = sqrt(stat.var);
  if (stat.var != 0.0) {
    stat.skew /= (double)(n * stat.var * stat.sdev);
    stat.curt = stat.curt / (double)(n * stat.var * stat.var) - 3.0;
  } else
    std::cerr << "@moment, No skew/kurtosis when variance = 0 (in moment)\n";
}


void gray_level_analysis() {
  std::cout << "GRAY LEVEL ANALYSIS" << std::endl;
  read_image(im_index_ref, iref, true); // Read reference image

  std::vector<std::ofstream> file(num_grains);
  for (int igrain = 0; igrain < num_grains; igrain++) {
    char fname[256];
    sprintf(fname, "gray_level%d.txt", igrain);
    file[igrain].open(fname);
  }

  for (int num_image = ibeg; num_image <= iend; num_image += iinc) {
    read_image(im_index_current, num_image);

    for (int igrain = 0; igrain < num_grains; igrain++) {

      std::vector<double> val;
      double min = 65535.0 + 1.0;
      double max = 0.0;
      int xpixel, ypixel;
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel = grain[igrain].pattern[i].dx + grain[igrain].refcoord_xpix;
        ypixel = grain[igrain].pattern[i].dy + grain[igrain].refcoord_ypix;
        double gray = (double)image[im_index_current][xpixel][ypixel];
        val.push_back(gray);
        if (min > gray)
          min = gray;
        if (max < gray)
          max = gray;
      }

      data_stat stat;
      moment(val, stat);

      file[igrain] << num_image << " " << stat.ave << " " << stat.adev << " " << stat.sdev << " " << stat.var << " "
                   << stat.skew << " " << stat.curt << " " << min << " " << max << '\n';
      // 1. image number
      // 2. average
      // 3. absolut deviation
      // 4. standard deviation
      // 5. variance
      // 6. skewness
      // 7. curtosis
      // 8. min
      // 9. max
    }
  }
}

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
  std::cout << "PATTERN QUALITY" << std::endl;
  read_image(0, iref, true); // Read image
  const int nbins = 100;
  int bin = 0;
  std::ofstream file("autocor.txt");
  std::vector<std::vector<double>> autoCor(num_grains);
  for (int i = 0; i < num_grains; i++)
    autoCor[i].resize(nbins);
  std::vector<std::vector<double>> population(num_grains);
  for (int i = 0; i < num_grains; i++)
    population[i].resize(nbins);
  double step = 0.05;
  double dist_max = sqrt(2.0 * 10.0 * 10.0);

  for (int igrain = 0; igrain < num_grains; igrain++) {

    std::cout << "Point " << igrain << std::endl;

    int x0 = grain[igrain].refcoord_xpix, y0 = grain[igrain].refcoord_ypix; // origin for auto-correlations
    double xmin = x0 - 10.0, xmax = x0 + 10.0;
    double ymin = y0 - 10.0, ymax = y0 + 10.0;

    for (int i = 0; i < nbins; i++) {
      autoCor[igrain][i] = 0.0;
      population[igrain][i] = 0;
    }

    std::vector<double> x1;
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

  std::ofstream fileMeanDev("autocorMeanDev.txt");
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
    fileMeanDev << dist << " " << Mean << " " << sqrt(Dev) << std::endl;
  }
}

/************************************************************************************************/
/*                                     POST PROCESSING                                          */
/************************************************************************************************/

// For post-processing specific to 1g2e see the fortran code in '1g2e_tools/postprocessDICfiles'
// Otherwise, for more generic processings, these procedures can be used.

struct sensor_data {
  double time; // seconds

  double L1; // cm
  double L2;
  double L3;
  double L4;

  double gamma; // degrees

  double sign1; // kPa
  double tau1;
  double sign2;
  double tau2;

  double fx20; // kN
  double fy20;
  double fx21;
  double fy21;
  double fx22;
  double fy22;
};

struct processed_point {
  double x, y, radius, rot;
  bool undistorded;
  bool scaled;
  processed_point() : x(0.0), y(0.0), radius(0.0), rot(0.0), undistorded(false), scaled(false) {}
};

struct process_data {
  char imFile_FROM[256];
  char imFile_TO[256];
  double scaleFactor;
  double Gxx, Gxy, Gyx, Gyy;
  double Exx, Exy, Eyx, Eyy;

  std::vector<Grain> gr;
  std::vector<sensor_data> sensor;

  std::vector<processed_point> pt;      // The points corresponding to sample elements
  std::vector<processed_point> corners; // 4 corners of the loading parallepiped
  std::vector<processed_point> fix_pt0; // 2 first fixe points on the first dic_out file
  std::vector<processed_point> fix_pt;  // All fixe points

  process_data() : scaleFactor(1.0) {}
};

void read_sensor_data(const char *name, process_data &data) {
  data.sensor.clear();
  std::ifstream file(name);
  // ignore 7 lines
  char ignoredLine[256];
  for (int i = 1; i <= 7; i++)
    file.getline(ignoredLine, 256);
  sensor_data S;

  while (file >> S.time >> S.L1 >> S.L2 >> S.L3 >> S.L4 >> S.gamma >> S.sign1 >> S.tau1 >> S.sign2 >> S.tau2 >>
         S.fx20 >> S.fy20 >> S.fx21 >> S.fy21 >> S.fx22 >> S.fy22) {
    data.sensor.push_back(S);
  }
}

//
void read(const char *name, process_data &data) {
  data.gr.clear();
  std::ifstream file(name);
  file >> num_grains;
  data.gr.resize(num_grains);
  for (int i = 0; i < num_grains; i++) {
    file >> data.gr[i].refcoord_xpix >> data.gr[i].refcoord_ypix >> data.gr[i].refrot >> data.gr[i].radius_pix >>
        data.gr[i].dx >> data.gr[i].dy >> data.gr[i].drot >> data.gr[i].upix >> data.gr[i].vpix >>
        data.gr[i].rot_inc >> data.gr[i].NCC >> data.gr[i].NCC_rescue >> data.gr[i].NCC_subpix;
  }
  if (file.good())
    file >> data.imFile_FROM >> data.imFile_TO;

  data.pt.clear();
  data.corners.clear();
  data.fix_pt.clear();

  for (size_t i = 0; i < data.gr.size(); i++) {
    processed_point P;
    P.x = data.gr[i].refcoord_xpix + data.gr[i].dx;
    P.y = data.gr[i].refcoord_ypix + data.gr[i].dy;
    P.radius = data.gr[i].radius_pix;
    P.rot = data.gr[i].refrot + data.gr[i].drot;

    if (P.radius == radius_corners)
      data.corners.push_back(P);
    else if (P.radius == radius_fixedPoints)
      data.fix_pt.push_back(P);
    else
      data.pt.push_back(P);
  }

  if (post_undistor) {
    double x, y, xu, yu;
    for (size_t i = 0; i < data.pt.size(); i++) {
      x = data.pt[i].x;
      y = data.pt[i].x;
      undistor(&disto_parameters[0], x, y, xu, yu);
      data.pt[i].undistorded = true;
      data.pt[i].x = xu;
      data.pt[i].y = yu;
    }
    for (size_t i = 0; i < data.corners.size(); i++) {
      x = data.corners[i].x;
      y = data.corners[i].x;
      undistor(&disto_parameters[0], x, y, xu, yu);
      data.corners[i].undistorded = true;
      data.corners[i].x = xu;
      data.corners[i].y = yu;
    }
    for (size_t i = 0; i < data.fix_pt.size(); i++) {
      x = data.fix_pt[i].x;
      y = data.fix_pt[i].x;
      undistor(&disto_parameters[0], x, y, xu, yu);
      data.fix_pt[i].undistorded = true;
      data.fix_pt[i].x = xu;
      data.fix_pt[i].y = yu;
    }
  }
}

// This function save the configuration after being post-processed
// in a file named CONFXXXX
void save_conf(int num, process_data &data) {
  char name[256];
  sprintf(name, "CONF%04d", num);
  std::ofstream conf(name);

  conf << data.pt.size() << std::endl;
  conf << std::scientific << std::setprecision(5);
  for (size_t i = 0; i < data.pt.size(); i++) {
    conf << data.pt[i].x << ' ' << data.pt[i].y << ' ' << data.pt[i].rot << ' ' << data.pt[i].radius << std::endl;
  }

  conf << "#CORNERS  The four corners" << std::endl;
  conf << data.corners.size() << std::endl;
  for (size_t i = 0; i < data.corners.size(); i++) {
    conf << data.corners[i].x << ' ' << data.corners[i].y << ' ' << data.corners[i].rot << std::endl;
  }

  conf << "#FIX_PTS  The fix points" << std::endl;
  conf << data.fix_pt.size() << std::endl;
  for (size_t i = 0; i < data.fix_pt.size(); i++) {
    conf << data.fix_pt[i].x << ' ' << data.fix_pt[i].y << ' ' << data.fix_pt[i].rot << std::endl;
  }

  conf << "#SCALES_  The scale factor expressed in length per pixel (e.g. m/pix)" << std::endl;
  conf << data.scaleFactor << std::endl;

  conf << "#FILES__  The two files used to correlate" << std::endl;
  conf << data.imFile_FROM << ' ' << data.imFile_TO << std::endl;
}

void post_process() {
  std::cout << "POST PROCESS" << std::endl;
  // dans le fichier de commande
  // ibeg <value>
  // iend <value>
  // iinc <value>
  // scaling_distance <value>
  // post_undistor <1|0>
  // post_rescale  <1|0>
  // post_sync <1|0>
  // post_rotation <1|0>
  // radius_corners <value>
  // radius_fixedPoints <value>

  if (post_sync) {
    // read sensor data from 1g2e
  }

  process_data data;

  // First file
  char name[256];
  sprintf(name, "dic_out_%d.txt", ibeg);
  read(name, data);

  if (post_rescale) {
    if (data.fix_pt.size() < 2) {
      std::cerr << "Cannot rescale, at least 2 fixed points need to be set" << std::endl;
      exit(0);
    }
    if (scaling_distance == 0.0) {
      std::cerr << "Cannot rescale, scaling_distance need to be set" << std::endl;
      exit(0);
    }
    data.fix_pt0.clear();
    data.fix_pt0.push_back(data.fix_pt[0]);
    data.fix_pt0.push_back(data.fix_pt[1]);
    double dx = data.fix_pt0[1].x - data.fix_pt0[0].x;
    double dy = data.fix_pt0[1].y - data.fix_pt0[0].y;
    double pix_length = sqrt(dx * dx + dy * dy);
    data.scaleFactor = scaling_distance / pix_length; // unit = length / pix
  }

  double transfo[3];
  transfo[0] = transfo[1] = transfo[2] = 0.0;

  int conf_num = 0;
  for (int num = ibeg; num <= iend; num += iinc) {
    char name[256];
    sprintf(name, "dic_out_%d.txt", num);
    read(name, data);

    if (post_rescale) {
      transfo[0] = data.fix_pt[0].x - data.fix_pt0[0].x;
      transfo[1] = data.fix_pt[0].y - data.fix_pt0[0].y;
      if (post_rotation) {
        double vx0 = data.fix_pt0[1].x - data.fix_pt0[0].x;
        double vy0 = data.fix_pt0[1].y - data.fix_pt0[0].y;
        double vx1 = data.fix_pt[1].x - data.fix_pt[0].x;
        double vy1 = data.fix_pt[1].y - data.fix_pt[0].y;
        double len0 = sqrt(vx0 * vx0 + vy0 * vy0);
        double len1 = sqrt(vx1 * vx1 + vy1 * vy1);
        double prod = vx0 * vy1 - vy0 * vx1;
        transfo[2] = asin(prod / (len0 * len1));
      } else {
        transfo[2] = 0.0;
      }

      double cosa = cos(-transfo[2]);
      double sina = sin(-transfo[2]);
      double xk, yk;

      // Global translation + rotation
      for (size_t i = 0; i < data.pt.size(); i++) {
        xk = data.pt[i].x - transfo[0];
        yk = transfo[1] - data.pt[i].y; // inversion of y-axis
        data.pt[i].x = xk * cosa - yk * sina;
        data.pt[i].y = xk * sina + yk * cosa;
      }
      for (size_t i = 0; i < data.fix_pt.size(); i++) {
        xk = data.fix_pt[i].x - transfo[0];
        yk = transfo[1] - data.fix_pt[i].y; // inversion of y-axis
        data.fix_pt[i].x = xk * cosa - yk * sina;
        data.fix_pt[i].y = xk * sina + yk * cosa;
      }
      for (size_t i = 0; i < data.fix_pt0.size(); i++) {
        xk = data.fix_pt0[i].x - transfo[0];
        yk = transfo[1] - data.fix_pt0[i].y; // inversion of y-axis
        data.fix_pt0[i].x = xk * cosa - yk * sina;
        data.fix_pt0[i].y = xk * sina + yk * cosa;
      }
      for (size_t i = 0; i < data.corners.size(); i++) {
        xk = data.corners[i].x - transfo[0];
        yk = transfo[1] - data.corners[i].y; // inversion of y-axis
        data.corners[i].x = xk * cosa - yk * sina;
        data.corners[i].y = xk * sina + yk * cosa;
      }

      // Rescaling
      for (size_t i = 0; i < data.pt.size(); i++) {
        data.pt[i].x *= data.scaleFactor;
        data.pt[i].y *= data.scaleFactor;
      }
      for (size_t i = 0; i < data.fix_pt.size(); i++) {
        data.fix_pt[i].x *= data.scaleFactor;
        data.fix_pt[i].y *= data.scaleFactor;
      }
      for (size_t i = 0; i < data.fix_pt0.size(); i++) {
        data.fix_pt0[i].x *= data.scaleFactor;
        data.fix_pt0[i].y *= data.scaleFactor;
      }
      for (size_t i = 0; i < data.corners.size(); i++) {
        data.corners[i].x *= data.scaleFactor;
        data.corners[i].y *= data.scaleFactor;
      }
    }

    save_conf(conf_num++, data);
  } // for loop over files
}

/************************************************************************************************/
/*                 DETERMINATION OF PARTICLE CENTERS WITH SUBPIXEL ACCURACY                     */
/************************************************************************************************/

double center_to_minimize_Circ(std::vector<double> &X) {
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
  std::cout << "FIND SUBPIXEL CENTERS" << std::endl;
  read_image(0, iref, true); // Read image

  // 3 profiles to identify the threshold value of gray
  std::ofstream prof("profile.txt");
  int yline = (int)(dimy / 4);
  for (int i = 0; i < dimx; i++) {
    prof << i << " " << image[0][i][yline] << std::endl;
  }
  prof << std::endl;
  yline = (int)(dimy / 2);
  for (int i = 0; i < dimx; i++) {
    prof << i << " " << image[0][i][yline] << std::endl;
  }
  prof << std::endl;
  yline = (int)(3 * dimy / 4);
  for (int i = 0; i < dimx; i++) {
    prof << i << " " << image[0][i][yline] << std::endl;
  }

  double x, y;
  double x0, y0;
  double x1, y1;
  double x2, y2;
  double x3, y3;
  double gray0, gray1, gray2, gray3;
  double gray, previous_gray;

  std::ofstream file("recentered.data");
  std::ofstream logFile("recentered.log.data");
  file << num_grains << std::endl << std::endl;
  logFile << num_grains << std::endl << std::endl;
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

    Powell<double(std::vector<double> &)> powell(center_to_minimize_Circ, 1e-8);
    center_parameters = powell.minimize(center_parameters, center_parameters_perturb);

    std::cout << "grain number " << igrain << std::endl;
    std::cout << "x offset = " << center_parameters[0] << std::endl;
    std::cout << "y offset = " << center_parameters[1] << std::endl;
    std::cout << "radius   = " << center_parameters[2] << std::endl;
    std::cout << std::endl;

    logFile << px0 << " " << py0 << " " << r0 << " " << center_parameters[0] << " " << center_parameters[1] << " "
            << center_parameters[2] << std::endl;
    file << px0 + center_parameters[0] << " " << py0 + center_parameters[1] << " " << grain[igrain].refrot << " "
         << center_parameters[2] << std::endl;
  } // for loop igrain
}

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
  std::cout << "def = (" << defx << ", " << defy << "),"
            << " trans = (" << transy_pix << ", " << transy_pix << "), rot = " << rot_deg << " deg. -> " << fname
            << std::endl;
  synth_im.writeTiff(fname, 0.0);
}

/************************************************************************************************/
/*                                        VISU PROCESS                                          */
/************************************************************************************************/

struct triangle_def {
  int i0, i1, i2;
  double exx, eyy, exy;
  double detJ;
};

struct point_def {
  double exx, eyy, exy;
  double detJ;
  point_def() : exx(0.0), eyy(0.0), exy(0.0), detJ(0.0) {}
};

// explanation can be found here: http://mathforum.org/library/drmath/view/72047.html
inline bool findQuadCoeff(double Xvalues[], double Yvalues[], double &a, double &b, double &c, int NumPoints) {
  double S00 = NumPoints;
  double S40 = 0.0, S10 = 0.0, S20 = 0.0, S30 = 0.0, S01 = 0.0, S11 = 0.0, S21 = 0.0;
  double SqrX;
  for (int i = 0; i < NumPoints; i++) {
    SqrX = Xvalues[i] * Xvalues[i];
    S40 += SqrX * SqrX;
    S30 += SqrX * Xvalues[i];
    S20 += SqrX;
    S10 += Xvalues[i];

    S01 += Yvalues[i];
    S11 += Xvalues[i] * Yvalues[i];
    S21 += SqrX * Yvalues[i];
  }

  double SqrS10 = S10 * S10;
  double SqrS20 = S20 * S20;
  double SqrS30 = S30 * S30;
  double denom = S00 * S20 * S40 - SqrS10 * S40 - S00 * SqrS30 + 2 * S10 * S20 * S30 - S20 * S20 * S20;
  if (denom == 0.0) {
    a = b = c = 0.0;
    return false;
  }
  double invDenom = 1.0 / denom;

  a = (S01 * S10 * S30 - S11 * S00 * S30 - S01 * SqrS20 + S11 * S10 * S20 + S21 * S00 * S20 - S21 * SqrS10) *
      invDenom;
  b = (S11 * S00 * S40 - S01 * S10 * S40 + S01 * S20 * S30 - S21 * S00 * S30 - S11 * SqrS20 + S21 * S10 * S20) *
      invDenom;
  c = (S01 * S20 * S40 - S11 * S10 * S40 - S01 * SqrS30 + S11 * S20 * S30 + S21 * S10 * S30 - S21 * SqrS20) *
      invDenom;

  return true;
}

void visu_process() {
  std::cout << "VISU PROCESS" << std::endl;

  int num_image;
  char fname[256];
  std::vector<Grain> grain_ref;

  sprintf(fname, dic_name, iref);
  fprintf(stdout, "\n____ ref DIC file %s\n", fname);
  read_grains(fname, true);
  for (size_t i = 0; i < grain.size(); i++)
    grain_ref.push_back(grain[i]);

  for (num_image = ibeg; num_image <= iend; num_image += iinc) {
    fprintf(stdout, "\n____ Image %d\n", num_image);
    read_image(0, num_image, true);

    sprintf(fname, dic_name, num_image);
    fprintf(stdout, "\n____ from DIC file %s\n", fname);
    read_grains(fname, true);

    thumbnail thumb(image, dimx, dimy, 0, image_div);
    ColorTable ct;
    ct.setMinMax(colorMin, colorMax);
    if (num_image == ibeg) {
      ct.savePpm("ColorMap.ppm");
    }
    colorRGBA cRGBA;

    if (visu_mode == "continuum") { /// ============================================================================
      std::cout << "CONTINUUM" << std::endl;
      std::vector<Point2D> points;
      for (int i = 0; i < num_grains; ++i) {
        double x_centre = grain[i].refcoord_xpix + grain[i].dx;
        double y_centre = grain[i].refcoord_ypix + grain[i].dy;

        points.push_back(Point2D(nearest(x_centre), nearest(y_centre)));
      }

      // MESH
      std::vector<triangle_def> triangles;
      if (!mesh.empty()) {
        triangle_def TD;
        for (size_t tr = 0; tr < mesh.size(); tr++) {
          TD.i0 = mesh[tr].i0;
          TD.i1 = mesh[tr].i1;
          TD.i2 = mesh[tr].i2;
          triangles.push_back(TD);
        }
      } else {
        Delaunay mygrid(points);
        triangle_def TD;
        for (size_t tr = 0; tr < mygrid.thelist.size(); tr++) {
          if (mygrid.thelist[tr].stat <= 0)
            continue;

          TD.i0 = mygrid.thelist[tr].p[0];
          TD.i1 = mygrid.thelist[tr].p[1];
          TD.i2 = mygrid.thelist[tr].p[2];
          triangles.push_back(TD);
        }
      }

      // Données aux noeuds
      std::vector<double> data;

      if (visu_output.substr(0, 6) == "strain") {
        // calculer les def par triangle ici
        double sumexx = 0.0;
        double sumeyy = 0.0;
        double sumexy = 0.0;
        double sumdetJ = 0.0;

        size_t i0, i1, i2;
        for (size_t tr = 0; tr < triangles.size(); tr++) {
          i0 = triangles[tr].i0;
          i1 = triangles[tr].i1;
          i2 = triangles[tr].i2;

          // valeurs
          double x1 = grain_ref[i0].refcoord_xpix + grain_ref[i0].dx;
          double y1 = grain_ref[i0].refcoord_ypix + grain_ref[i0].dy;
          double x2 = grain_ref[i1].refcoord_xpix + grain_ref[i1].dx;
          double y2 = grain_ref[i1].refcoord_ypix + grain_ref[i1].dy;
          double x3 = grain_ref[i2].refcoord_xpix + grain_ref[i2].dx;
          double y3 = grain_ref[i2].refcoord_ypix + grain_ref[i2].dy;

          double ux1 = grain[i0].refcoord_xpix + grain[i0].dx - x1;
          double uy1 = grain[i0].refcoord_ypix + grain[i0].dy - y1;
          double ux2 = grain[i1].refcoord_xpix + grain[i1].dx - x2;
          double uy2 = grain[i1].refcoord_ypix + grain[i1].dy - y2;
          double ux3 = grain[i2].refcoord_xpix + grain[i2].dx - x3;
          double uy3 = grain[i2].refcoord_ypix + grain[i2].dy - y3;

          // calcul de def
          double x13 = x1 - x3;
          double y13 = y1 - y3;
          double x23 = x2 - x3;
          double y23 = y2 - y3;

          triangles[tr].detJ = fabs(x13 * y23 - y13 * x23); // det J = 2 * Surf
          sumdetJ += triangles[tr].detJ;

          double x32 = x3 - x2;
          double x21 = x2 - x1;
          double y31 = y3 - y1;
          double y12 = y1 - y2;

          // This is not exactly strain but (2 * Surf * strain)
          triangles[tr].exx = y23 * ux1 + y31 * ux2 + y12 * ux3;
          triangles[tr].eyy = x32 * uy1 + x13 * uy2 + x21 * uy3;
          triangles[tr].exy = x32 * ux1 + y23 * uy1 + x13 * ux2 + y31 * uy2 + x21 * ux3 + y12 * uy3;

          sumexx += triangles[tr].exx;
          sumeyy += triangles[tr].eyy;
          sumexy += triangles[tr].exy;
        }

        double exx_mean = sumexx;
        double eyy_mean = sumeyy;
        double exy_mean = sumexy;
        if (sumdetJ > 0.0) {
          exx_mean /= sumdetJ;
          eyy_mean /= sumdetJ;
          exy_mean /= sumdetJ;
        }

        std::cout << "exx_mean = " << exx_mean << std::endl;
        std::cout << "eyy_mean = " << eyy_mean << std::endl;
        std::cout << "exy_mean = " << exy_mean << std::endl;

        std::vector<point_def> point_defs;
        point_def PD;
        for (size_t i = 0; i < grain.size(); ++i)
          point_defs.push_back(PD);
        for (size_t t = 0; t < triangles.size(); t++) {
          point_defs[triangles[t].i0].exx += triangles[t].exx;
          point_defs[triangles[t].i0].eyy += triangles[t].eyy;
          point_defs[triangles[t].i0].exy += triangles[t].exy;
          point_defs[triangles[t].i0].detJ += triangles[t].detJ;

          point_defs[triangles[t].i1].exx += triangles[t].exx;
          point_defs[triangles[t].i1].eyy += triangles[t].eyy;
          point_defs[triangles[t].i1].exy += triangles[t].exy;
          point_defs[triangles[t].i1].detJ += triangles[t].detJ;

          point_defs[triangles[t].i2].exx += triangles[t].exx;
          point_defs[triangles[t].i2].eyy += triangles[t].eyy;
          point_defs[triangles[t].i2].exy += triangles[t].exy;
          point_defs[triangles[t].i2].detJ += triangles[t].detJ;
        }

        for (size_t i = 0; i < point_defs.size(); i++) {
          if (point_defs[i].detJ > 0.0) {
            point_defs[i].exx /= point_defs[i].detJ;
            point_defs[i].eyy /= point_defs[i].detJ;
            point_defs[i].exy /= point_defs[i].detJ;
          }

          if (visu_output == "strainxx") {
            data.push_back(point_defs[i].exx);
          } else if (visu_output == "strainyy") {
            data.push_back(point_defs[i].eyy);
          } else if (visu_output == "strainxy") {
            data.push_back(point_defs[i].exy);
          } else {
            data.push_back(0.0);
          }
        }
      } else {
        for (int i = 0; i < num_grains; ++i) {
          if (visu_output == "NCC") {
            data.push_back(grain[i].NCC);
          } else if (visu_output == "NCC_rescue") {
            data.push_back(grain[i].NCC_rescue);
          } else if (visu_output == "NCC_subpix") {
            data.push_back(grain[i].NCC_subpix);
          } else if (visu_output == "disp") {
            double dx = (grain[i].refcoord_xpix + grain[i].dx) - (grain_ref[i].refcoord_xpix + grain_ref[i].dx);
            double dy = (grain[i].refcoord_ypix + grain[i].dy) - (grain_ref[i].refcoord_ypix + grain_ref[i].dy);
            data.push_back(sqrt(dx * dx + dy * dy));
          } else if (visu_output == "dispx") {
            double dx = (grain[i].refcoord_xpix + grain[i].dx) - (grain_ref[i].refcoord_xpix + grain_ref[i].dx);
            data.push_back(dx);
          } else if (visu_output == "dispy") {
            double dy = (grain[i].refcoord_ypix + grain[i].dy) - (grain_ref[i].refcoord_ypix + grain_ref[i].dy);
            data.push_back(dy);
          } else {
            data.push_back(0.0);
          }
        }
      }

      if (visu_autoscale == 1) {
        colorMin = data[0];
        colorMax = data[0];
        for (int i = 1; i < num_grains; ++i) {
          if (data[i] < colorMin)
            colorMin = data[i];
          if (data[i] > colorMax)
            colorMax = data[i];
        }
        ct.setMinMax(colorMin, colorMax);
      }

      size_t i0, i1, i2;

      RGB col0, col1, col2;
      for (size_t tr = 0; tr < triangles.size(); tr++) {
        i0 = triangles[tr].i0;
        i1 = triangles[tr].i1;
        i2 = triangles[tr].i2;

        ct.getRGB(data[i0], &cRGBA);
        col0.r = cRGBA.r;
        col0.g = cRGBA.g;
        col0.b = cRGBA.b;
        ct.getRGB(data[i1], &cRGBA);
        col1.r = cRGBA.r;
        col1.g = cRGBA.g;
        col1.b = cRGBA.b;
        ct.getRGB(data[i2], &cRGBA);
        col2.r = cRGBA.r;
        col2.g = cRGBA.g;
        col2.b = cRGBA.b;

        double sc = 1.0 / image_div;
        thumb.draw_triangle((int)(points[i0].x * sc), (int)(points[i0].y * sc), (int)(points[i1].x * sc),
                            (int)(points[i1].y * sc), (int)(points[i2].x * sc), (int)(points[i2].y * sc), col0, col1,
                            col2);

        if (0) {
          RGB yellow(255, 255, 0);
          thumb.draw_line((int)(points[i0].x * sc), (int)(points[i0].y * sc), (int)(points[i1].x * sc),
                          (int)(points[i1].y * sc), yellow);
          thumb.draw_line((int)(points[i1].x * sc), (int)(points[i1].y * sc), (int)(points[i2].x * sc),
                          (int)(points[i2].y * sc), yellow);
          thumb.draw_line((int)(points[i2].x * sc), (int)(points[i2].y * sc), (int)(points[i0].x * sc),
                          (int)(points[i0].y * sc), yellow);
        }
      }
    } else if (visu_mode == "grid") { /// ========================================================================

      int ny_grid_visu = num_grains / nx_grid_visu;
      int npts = 2 * smoothDegreeF + 1;

      // Some checks
      if (num_grains % nx_grid_visu != 0) {
        std::cerr << "nx_grid_visu not compatible with num_grains" << std::endl;
        return;
      }
      if (nx_grid_visu < npts) {
        std::cerr << "nx_grid_visu too small" << std::endl;
        return;
      }
      if (ny_grid_visu < npts) {
        std::cerr << "ny_grid_visu too small" << std::endl;
        return;
      }

      std::vector<double> Fxx(num_grains), Fxy(num_grains), Fyx(num_grains), Fyy(num_grains);
      std::vector<double> Exx(num_grains), Exy(num_grains), Eyx(num_grains), Eyy(num_grains);
      std::vector<double> XData(npts), YData(npts), UxData(npts), UyData(npts);
      std::vector<double> data;

      if (visu_output.substr(0, 6) == "strain") {

        // Compute transformation gradient matrix
        for (int j = 0; j < ny_grid_visu; j++) {
          for (int i = 0; i < nx_grid_visu; i++) {
            int ig = j * nx_grid_visu + i;

            int i0 = i - smoothDegreeF;
            if (i0 < 0)
              i0 = 0;
            if (i0 > nx_grid_visu - npts)
              i0 = nx_grid_visu - npts;
            double x0 = grain_ref[ig].refcoord_xpix + grain_ref[ig].dx;
            for (int n = i0; n < i0 + npts; n++) {
              int ii = j * nx_grid_visu + n;
              XData[n - i0] = (grain_ref[ii].refcoord_xpix + grain_ref[ii].dx) - x0;
              UxData[n - i0] =
                  (grain[ii].refcoord_xpix + grain[ii].dx) - (grain_ref[ii].refcoord_xpix + grain_ref[ii].dx);
            }

            int j0 = j - smoothDegreeF;
            if (j0 < 0)
              j0 = 0;
            if (j0 > ny_grid_visu - npts)
              j0 = ny_grid_visu - npts;
            double y0 = grain_ref[ig].refcoord_ypix + grain_ref[ig].dy;
            for (int n = j0; n < j0 + npts; n++) {
              int ii = n * nx_grid_visu + i;
              YData[n - j0] = (grain_ref[ii].refcoord_ypix + grain_ref[ii].dy) - y0;
              UyData[n - j0] =
                  grain[ii].refcoord_ypix + grain[ii].dy - (grain_ref[ii].refcoord_ypix + grain_ref[ii].dy);
            }

            double A, B, C;
            bool success;
            success = findQuadCoeff(&XData[0], &UxData[0], A, B, C, npts);
            if (success)
              Fxx[ig] = B + 1.0;
            else
              Fxx[ig] = 0.0;
            success = findQuadCoeff(&YData[0], &UyData[0], A, B, C, npts);
            if (success)
              Fyy[ig] = B + 1.0;
            else
              Fyy[ig] = 0.0;
            success = findQuadCoeff(&YData[0], &UxData[0], A, B, C, npts);
            if (success)
              Fxy[ig] = B;
            else
              Fxy[ig] = 0.0;
            success = findQuadCoeff(&XData[0], &UyData[0], A, B, C, npts);
            if (success)
              Fyx[ig] = B;
            else
              Fyx[ig] = 0.0;
          }
        }

        // Compute the Green-Lagrange Strain tensor E = (FFt - I)/2
        for (int i = 0; i < num_grains; ++i) {
          Exx[i] = 0.5 * (Fxx[i] * Fxx[i] + Fyx[i] * Fyx[i] - 1.0);
          Exy[i] = 0.5 * (Fxx[i] * Fxy[i] + Fyx[i] * Fyy[i]);
          Eyx[i] = 0.5 * (Fxy[i] * Fxx[i] + Fyy[i] * Fyx[i]);
          Eyy[i] = 0.5 * (Fxy[i] * Fxy[i] + Fyy[i] * Fyy[i] - 1.0);
        }

        // Compute data asked by the user
        for (int i = 0; i < num_grains; ++i) {
          if (visu_output == "strainxx") {
            data.push_back(Exx[i]);
          } else if (visu_output == "strainyy") {
            data.push_back(Eyy[i]);
          } else if (visu_output == "strainxy") {
            data.push_back(Exy[i]);
          } else {
            data.push_back(0.0);
          }
        }
      } // End strain computation
      else { // no strain computation is required

        // Compute data asked by the user
        for (int i = 0; i < num_grains; ++i) {
          if (visu_output == "NCC") {
            data.push_back(grain[i].NCC);
          } else if (visu_output == "NCC_rescue") {
            data.push_back(grain[i].NCC_rescue);
          } else if (visu_output == "NCC_subpix") {
            data.push_back(grain[i].NCC_subpix);
          } else if (visu_output == "dispx") {
            double dx = (grain[i].refcoord_xpix + grain[i].dx) - (grain_ref[i].refcoord_xpix + grain_ref[i].dx);
            data.push_back(dx);
          } else if (visu_output == "dispy") {
            double dy = (grain[i].refcoord_ypix + grain[i].dy) - (grain_ref[i].refcoord_ypix + grain_ref[i].dy);
            data.push_back(dy);
          } else if (visu_output == "disp") {
            double dx = (grain[i].refcoord_xpix + grain[i].dx) - (grain_ref[i].refcoord_xpix + grain_ref[i].dx);
            double dy = (grain[i].refcoord_ypix + grain[i].dy) - (grain_ref[i].refcoord_ypix + grain_ref[i].dy);
            data.push_back(sqrt(dx * dx + dy * dy));
          } else {
            data.push_back(0.0);
          }
        }
      }

      // Set the color range
      if (visu_autoscale == 1) {
        std::cout << "Rescaling the colormap" << std::endl;
        colorMin = data[0];
        colorMax = data[0];
        for (int i = 1; i < num_grains; ++i) {
          if (data[i] < colorMin)
            colorMin = data[i];
          if (data[i] > colorMax)
            colorMax = data[i];
        }
        ct.setMinMax(colorMin, colorMax);
        TRACKER_LOG("ColorScale / " << visu_output << " / im = " << num_image << " / min = " << colorMin
                                    << " / max = " << colorMax << std::endl);
      }

      int p0x, p1x, p2x, p3x;
      int p0y, p1y, p2y, p3y;
      int ig0, ig1, ig2, ig3;

      RGB col0, col1, col2, col3;

      for (int ey = 0; ey < ny_grid_visu - 1; ey++) {
        for (int ex = 0; ex < nx_grid_visu - 1; ex++) {
          p0x = ex;
          p0y = ey;
          p1x = ex + 1;
          p1y = ey + 1;
          p2x = ex;
          p2y = ey + 1;
          p3x = ex + 1;
          p3y = ey;

          ig0 = p0y * nx_grid_visu + p0x;
          ig1 = p1y * nx_grid_visu + p1x;
          ig2 = p2y * nx_grid_visu + p2x;
          ig3 = p3y * nx_grid_visu + p3x;

          ct.getRGB(data[ig0], &cRGBA);
          col0.r = cRGBA.r;
          col0.g = cRGBA.g;
          col0.b = cRGBA.b;
          ct.getRGB(data[ig1], &cRGBA);
          col1.r = cRGBA.r;
          col1.g = cRGBA.g;
          col1.b = cRGBA.b;
          ct.getRGB(data[ig2], &cRGBA);
          col2.r = cRGBA.r;
          col2.g = cRGBA.g;
          col2.b = cRGBA.b;
          ct.getRGB(data[ig3], &cRGBA);
          col3.r = cRGBA.r;
          col3.g = cRGBA.g;
          col3.b = cRGBA.b;

          double sc = 1.0 / (double)image_div;

          if (visu_draw_in_ref == 1) {
            thumb.draw_triangle((grain[ig0].refcoord_xpix) * sc, (grain[ig0].refcoord_ypix) * sc,
                                (grain[ig1].refcoord_xpix) * sc, (grain[ig1].refcoord_ypix) * sc,
                                (grain[ig2].refcoord_xpix) * sc, (grain[ig2].refcoord_ypix) * sc, col0, col1, col2);
            thumb.draw_triangle((grain[ig0].refcoord_xpix) * sc, (grain[ig0].refcoord_ypix) * sc,
                                (grain[ig3].refcoord_xpix) * sc, (grain[ig3].refcoord_ypix) * sc,
                                (grain[ig1].refcoord_xpix) * sc, (grain[ig1].refcoord_ypix) * sc, col0, col3, col1);
          } else {
            thumb.draw_triangle(
                (grain[ig0].refcoord_xpix + grain[ig0].dx) * sc, (grain[ig0].refcoord_ypix + grain[ig0].dy) * sc,
                (grain[ig1].refcoord_xpix + grain[ig1].dx) * sc, (grain[ig1].refcoord_ypix + grain[ig1].dy) * sc,
                (grain[ig2].refcoord_xpix + grain[ig2].dx) * sc, (grain[ig2].refcoord_ypix + grain[ig2].dy) * sc,
                col0, col1, col2);
            thumb.draw_triangle(
                (grain[ig0].refcoord_xpix + grain[ig0].dx) * sc, (grain[ig0].refcoord_ypix + grain[ig0].dy) * sc,
                (grain[ig3].refcoord_xpix + grain[ig3].dx) * sc, (grain[ig3].refcoord_ypix + grain[ig3].dy) * sc,
                (grain[ig1].refcoord_xpix + grain[ig1].dx) * sc, (grain[ig1].refcoord_ypix + grain[ig1].dy) * sc,
                col0, col3, col1);
          }
        }
      }
    } else if (visu_mode == "vectors") { /// ========================================================================

    } else if (visu_mode == "color_disks") { /// ====================================================================
      std::cout << "COLOR_DISKS" << std::endl;
      RGB col;
      std::vector<double> data;
      for (int i = 0; i < num_grains; ++i) {
        if (visu_output == "NCC") {
          data.push_back(grain[i].NCC);
        } else if (visu_output == "NCC_rescue") {
          data.push_back(grain[i].NCC_rescue);
        } else if (visu_output == "NCC_subpix") {
          data.push_back(grain[i].NCC_subpix);
        } else if (visu_output == "disp") {
          double dx = grain[i].refcoord_xpix + grain[i].dx - (grain_ref[i].refcoord_xpix + grain_ref[i].dx);
          double dy = grain[i].refcoord_ypix + grain[i].dy - (grain_ref[i].refcoord_ypix + grain_ref[i].dy);
          data.push_back(sqrt(dx * dx + dy * dy));
        } else {
          data.push_back(0.0);
        }
      }

      if (visu_autoscale) {
        std::cout << "Rescaling the colormap" << std::endl;
        colorMin = data[0];
        colorMax = data[0];
        for (int i = 1; i < num_grains; ++i) {
          if (data[i] < colorMin)
            colorMin = data[i];
          if (data[i] > colorMax)
            colorMax = data[i];
        }
        TRACKER_SHOW(colorMin);
        TRACKER_SHOW(colorMax);
        ct.setMinMax(colorMin, colorMax);
      }

      for (int i = 0; i < num_grains; ++i) {
        int x_centre = (int)(grain[i].refcoord_xpix + grain[i].dx);
        int y_centre = (int)(grain[i].refcoord_ypix + grain[i].dy);
        int rayon = (int)grain[i].radius_pix;

        x_centre /= image_div;
        y_centre /= image_div;
        rayon /= image_div;

        ct.getRGB(data[i], &cRGBA);
        col.r = cRGBA.r;
        col.g = cRGBA.g;
        col.b = cRGBA.b;
        thumb.draw_disk(x_centre, y_centre, rayon, col);
      }
    }

    sprintf(fname, "Visu_%s_%d.ppm", visu_output.c_str(), num_image);
    // thumb.writePpm(fname, visu_alpha);
    thumb.writeTiff(fname, visu_alpha);
  } // End loop over files (image and dic_out_x.txt)
}
