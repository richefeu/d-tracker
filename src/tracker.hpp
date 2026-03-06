// TRACKER
// December 2011, Janvier 2012
// Gael.Combe@ujf-grenoble.fr
// Vincent.Richefeu@ujf-grenoble.fr
// Lab 3SR, Grenoble

#include <cmath>
#include <cstdio>
#include <cstdlib>
#ifndef M_PI
#define M_PI 3.1415926535
#endif
#include <algorithm>
#include <cstdint> // for uint16_t
#include <cstring> // for strcmp
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

//#include <fftw3.h>

#define MAX_NUMBER_THREADS 48

#if defined(_OPENMP)
#include <omp.h>
#endif

int progress;

#include <tiffio.h>

#ifndef CYGWIN
#include "libraw.h"
#endif

#include "ColorTable.cpp"

using namespace std;

#define SHOW(V) cout << #V " = " << V << flush << endl

ofstream g_logfile;
#define LOG(V)                                                                                                       \
  cout << V;                                                                                                         \
  if (g_logfile)                                                                                                     \
  g_logfile << V

// const double convDepth = 1.0;

// All functions

// _______INLINED FONCTIONS________
// convert a double or float to the nearest integer
inline int nearest(double x) { return (int)floor(x + 0.5); }
inline int nearest(float x) { return (int)floor(x + 0.5); }

// fonction max
inline double max(double a, double b) { return ((a > b) ? a : b); }
inline double max(float a, float b) { return ((a > b) ? a : b); }

void header();
int init(int argc, char *argv[]);
void rotate_pixel_pattern(int igrain, int i, double c, double s, int *xpixel, int *ypixel);

double get_time();

void do_precomputations();
void follow_pattern_pixel(int igrain);
void follow_pattern_rescue_pixel(int igrain);
void follow_pattern_super_rescue_pixel(int igrain);

void follow_pattern_subpixel_xyR(int igrain);
double NCC_to_minimize_xyR(vector<double> &);

double Equiproj_Error_Computation(vector<double> &);
double disto_to_minimize_Equiproj(vector<double> &);
double disto_to_minimize_grid(vector<double> &);

double SolidMotionError_to_minimize(vector<double> &X);
void SolidMotionError(double &angle, double &xc, double &yc, double &xtrans, double &ytrans, double &mean,
                      double &stddev);

int read_data(const char *name); // To read the command file
void write_data(const char *name);
void read_image(int i, int num, bool = false);
void read_raw_image(int i, int num, bool = false);
void make_grid(double xmin, double xmax, double ymin, double ymax, int nx, int ny);
int read_grains(const char *name, bool restart);
int read_gmsh(const char *name);
int save_grains(int num);
int save_grains(const char *name, int num, bool simpleVersion = false);
void create_image(int num);
void create_image_Netbpm(int num);
void find_neighbours(int igrain);
void make_circ_pattern(int igrain, int radius_pix);
void make_ring_pattern(int igrain, int radius_IN_pix, int radius_OUT_pic);
void make_rect_pattern(int igrain, int half_width_pix, int half_height_pix);
void make_custom_pattern(const char *name);
void distor(double *X, const double xd, const double yd, double &xu, double &yu);
void diff_distor(double *X, const double xu, const double yu, double &dxd_dxu, double &dyd_dyu);
int undistor(double *X, const double xd, const double yd, double &xu, double &yu, int niterMax = 100,
             double tol = 1e-3, int nconv = 3);
void setSolidSolutionFromRef(double DX, double DY, double DROT);

void set_solidTranslation(double u, double v);
void mask_rect(int xmin, int xmax, int ymin, int ymax);

// Procedures
void particle_tracking();
void particle_tracking_v2_subpix();
void correction_distortion();
void correction_distortion_grid();
void find_subpixel_centers();
void pattern_quality();
void pattern_fft();
void post_process();

// Synthetic images
inline double findnoise2(double x, double y);
inline double interpolate(double a, double b, double x);
double noise(double x, double y);
void generate_synthetic_images(int w, int h, double zoom, int octaves, double p, double defx, double defy,
                               double transx_pix, double transy_pix, double rot_deg);

// ========================================================================================

struct relative_coord_type {
  int dx, dy; // Relative positions can be negative
};

struct subpix_coord_type {
  double x, y;
};

string procedure("particle_tracking"); // the default procedure is to track particles

int fake_undistor = 0;

// Distortion parameters (grid method)
string grid_image_name;
int nx_grid_disto;
int ny_grid_disto;

// Distortion parameters (displacement method)
vector<int> image_numbers_corrDisto;
vector<double> imposed_displ_x_pix, imposed_displ_y_pix;
vector<vector<double>> dx_corrDisto, dy_corrDisto, NCC_subpix_corrDisto;
vector<double> disto_parameters(8);
vector<double> disto_parameters_perturb(8);
// 0 -> xc
// 1 -> yc
// 2 -> K1
// 3 -> K2
// 4 -> K3
// 5 -> P1
// 6 -> P2
// 7 -> P3

// These vectors are used for the correction of distortions
vector<int> List_Grains_i;
vector<int> List_Grains_j;
double equiProj_dist_min = 0.0;
double equiProj_dist_max = 1000.0;

// subpixel positionning of rod centers
vector<double> xbound, ybound;
vector<double> center_parameters(3);
vector<double> center_parameters_perturb(3);
double subpix_center_dr0 = 10.0;
double subpix_center_xstep = 0.1;
double subpix_center_ystep = 0.1;
double subpix_center_rstep = 0.1;
int subpix_center_threshold = 10000;

// this allows to define special pattern according to a specified radius
// A negative radius means "all the grains"
double targetRadiusPattern = -1.0;

// Pattern quality
int halfPatternQual = 9;

// Post-Process
int post_undistor = 1;
int post_rescale = 1;
int post_sync = 0;
int post_rotation = 1;
char device_data_name[256];
double radius_corners = 1.0;
double radius_fixedPoints = 2.0;
double scaling_distance = 0.0;

int ShowSolidMotionError = 0;

// Visualization
float colorMin = 0.7;
float colorMax = 1.0;
int visu_autoscale = 1;
double visu_alpha = 0.3;
string visu_output = "NCC";
string visu_mode = "discrete";
int nx_grid_visu;
int smoothDegreeF = 1;
int visu_draw_in_ref = 0;

// Minimization (Subpixel)
#include "powell_nr3.hpp"
double subpix_tol = 1e-8;
double initial_direction_dx = 0.01;
double initial_direction_dy = 0.01;
double initial_direction_dr = 0.01;

double subpix_displacement_max = 3.0;
double overflow_stiffness = 0.0;

// Global variables

int RawImages = 0;

int DemosaicModel = 0;
int rescaleGrayLevels = 0;

struct imageDATA {
  float iso_speed;
  float shutter;
  float aperture;
  float focal_len;
  unsigned shot_order;
  std::string dateTime;
};

vector<vector<vector<uint16_t>>> image;   // image[0 or 1][dimx][dimy]
vector<imageDATA> imageData(2);           // imageData[0 or 1]
char image_name[256];                     // Image name in "printf" style. Example: "./quelquepart/toto%04d.tif"
char dic_name[256];                       // DIC-file name in "printf" style. Example: "./quelquepart/dic_out_%d.txt"
int iref, ibeg, iend, iinc, iraz, idelta; // Control of file numbers
bool require_precomputations = true;      // If true, remake the pre-computations
int dimx = 0, dimy = 0;                   // Image sizes
string pattern_command;
string grain_positions_command;

int im_index_ref = 0, im_index_current = 1; // Index of image origin and target in table image
int wanted_num_threads = 4;

int num_neighbour_max = 40;
double neighbour_dist_pix = 250.0; // max distance between points to find neighbours
int period_rebuild_neighbour_list = 50;

// Flags
int verbose_level = 0;
int rescue_level = 2; // in range [0 2].
                      // A negative value means no pixel-resolution tracking (go directly to subpixel optimization)
int subpixel = 1;
int rotations = 1;
int use_neighbour_list = 1;

// Image generation for check
int make_images = 0;     // create images (1) or not (0)
double image_size = 0.2; // in range [0 1] (DEPRECATED)
int image_div = 1;       // will replace image_size (TODO)
int draw_angle = 1;      //
int draw_disp = 1;       //
int draw_rescued = 1;    //

// ========================================================================================

// Process has done i out of n rounds,
// and we want a bar of width w and resolution r.
static inline void loadbar(size_t x, size_t n, size_t w = 50) {
  if ((x != n) && (x % (n / 100 + 1) != 0))
    return;

  float ratio = x / (float)n;
  size_t c = ratio * w;

  cerr << setw(3) << (size_t)(ratio * 100) << "% [";
  for (size_t x = 0; x < c; x++)
    cerr << "|";
  for (size_t x = c; x < w; x++)
    cerr << " ";
  cerr << "]\r" << flush;
}

std::string timestamp2string(time_t rawtime) {
  struct tm * timeinfo;
  timeinfo = localtime (&rawtime);
  char retChar[256];
  sprintf(retChar, "%s", asctime(timeinfo));
  return std::string(retChar);
}

#include "interpol.cpp"

struct grain_type_2D {
  // variables initiales (generalement ce sont les donnees d'entree)
  int refcoord_xpix; // coordonnée x de référence du point suivi (centre de masse)
  int refcoord_ypix; // coordonnée y de référence du point suivi (centre de masse)
  double refrot;     // angle initial (de référence) du grain
  double radius_pix; // Grain radius expressed in sub-pixels (input data)

  // Deplacements et rotation cumules par rapport aux variables initiales (par rapport à refxxx)
  double dx=0;
  double dy=0;
  double drot=0;

  // Backup des deplacements et rotations cumulés
  double dx_prev=0;
  double dy_prev=0;
  double drot_prev=0;

  // Kinematic data obtain on current DIC (num_image to num_image+iinc)
  double upix;    // x translation resulting from last DIC
  double vpix;    // y translation resulting from last DIC
  double rot_inc; // Increment of rotation resulting from last DIC

  // Coefficients de correlation
  double NCC;         // Best NCC
  double NCC_rescue;  // Best NCC after rescue (0 if no rescue has been performed)
  double NCC_subpix;  // Best NCC after subpix (0 if no subpix refinement has been performed)
  double mean0, C0C0; // Precomputed data for the computation of NCCs
  bool masked;        // Pour les point dans une zone à ne pas correlé (utilisation avec une grille)

  vector<relative_coord_type> pattern;          // A list of relative positions that represent a pattern
  vector<relative_coord_type> pattern0_rotated; // Rotated pattern for the reference position
  vector<relative_coord_type> pattern1_rotated; // Rotated pattern for the tested position

  vector<int> neighbour; // List of neighbouring grains
  int num_neighbour;     // Number of neighbour for igrain (according to neighbour_dist_pix) @fixme supprimer

  // Constructor
  grain_type_2D() {
    refcoord_xpix = 0;
    refcoord_ypix = 0;
    refrot = 0.0;
    radius_pix = 0.0;

    dx = dy = drot = 0.0;
    dx_prev = dy_prev = drot_prev = 0.0;
    upix = vpix = rot_inc = 0.0;
    NCC = NCC_rescue = NCC_subpix = 0.0;
    mean0 = C0C0 = 0.0;
    masked = false;
    num_neighbour = 0;
  }

  // Copy constructor
  grain_type_2D(const grain_type_2D &c) {
    refcoord_xpix = c.refcoord_xpix;
    refcoord_ypix = c.refcoord_ypix;
    refrot = c.refrot;
    radius_pix = c.radius_pix;

    dx = c.dx;
    dy = c.dy;
    drot = c.drot;

    dx_prev = c.dx_prev;
    dy_prev = c.dy_prev;
    drot_prev = c.drot_prev;

    upix = c.upix;
    vpix = c.vpix;
    rot_inc = c.rot_inc;
    NCC = c.NCC;
    NCC_rescue = c.NCC_rescue;
    NCC_subpix = c.NCC_subpix;

    mean0 = c.mean0;
    C0C0 = c.C0C0;
    masked = c.masked;
    num_neighbour = c.num_neighbour;

    // other data are not copied
  }

  // Destructor
  ~grain_type_2D() {}

  void backup() {
    dx_prev = dx;
    dy_prev = dy;
    drot_prev = drot;
  }

  void reset() { NCC = NCC_rescue = NCC_subpix = 0.0; }
};

struct triangle {
  int i0, i1, i2;
};

struct increasingRadius_order {
  inline bool operator()(const grain_type_2D &G1, const grain_type_2D &G2) { return (G1.radius_pix < G2.radius_pix); }
};

struct decreasingRadius_order {
  inline bool operator()(const grain_type_2D &G1, const grain_type_2D &G2) { return (G1.radius_pix > G2.radius_pix); }
};

struct increasingHeight_order {
  inline bool operator()(const grain_type_2D &G1, const grain_type_2D &G2) {
    return (G1.refcoord_ypix < G2.refcoord_ypix);
  }
};

vector<grain_type_2D> grain;
int num_grains = 0;

vector<triangle> mesh;

vector<int> to_be_rescued;
int num_to_be_rescued = 0;
double NCC_min = 0.7; // The NCC-value above which a first rescue is mandatory

vector<int> to_be_super_rescued;
int num_to_be_super_rescued = 0;
double NCC_min_super = 0.5; // The NCC-value above which a super_rescue is mandatory

// Research zone (for pixel-precision tracking)
struct search_zone_type {
  // Positive values
  int left;       // Towards decreasing x
  int right;      // Towards increasing x
  int up;         // Towards decreasing y
  int down;       // Towards increasing y
  int num_rot;    // Number of tested angles (above and after, total = 2*num_rot+1)
  double inc_rot; // Increment of rotation
};

search_zone_type search_zone;              // For non sub-pixel tracking
search_zone_type search_zone_rescue;       // For points with bad mathing
search_zone_type search_zone_super_rescue; // For points with very bad mathing

// A table that hold the grain number being processed for each thread.
// It is required in subpixel tracking.
// Do not hesitate to increase the maximum number of threads if necessary
int igrain_of_thread[MAX_NUMBER_THREADS];

#include "delaunay2D.hpp"
#include "image_io.cpp"
#include "thumbnail.hpp"
