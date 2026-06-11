/*
 * Tracker - Organisation du code
 *
 * Dans cette version, le code a été consolidé dans un seul fichier (tracker.cpp)
 * pour simplifier le développement.
 *
 * Choix techniques :
 * - Réduction du nombre de fichiers sources pour une maintenance plus aisée.
 * - Conservation de l'existant éprouvé, même si cela va à l'encontre de certaines
 *   "bonnes pratiques" de programmation.
 *
 * Compatibilité :
 * - Le code est globalement écrit en C++17, mais des éléments de style C++ ancien subsistent.
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint> // for uint16_t
#include <cstdio>
#include <cstdlib>
#include <cstring> // for strcmp
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#define MAX_NUMBER_THREADS 48

// interpol.hpp est header-only (methodes inline) et fournit les types des
// interpolateurs utilises par les variables globales declarees plus bas. Il peut
// donc etre inclus en toute securite par plusieurs unites de compilation.
#include "interpol.hpp"

// NOTE: Les includes "lourds" (libraw, tiffio, et les en-tetes toofus qui
// definissent des fonctions non-inline au niveau namespace) sont volontairement
// places dans tracker.cpp uniquement. Les inclure ici les ferait apparaitre dans
// chaque unite incluant ce header (run.cpp, tracker_gui.cpp) -> symboles dupliques
// a l'edition de liens avec libtracker.a.

#define TRACKER_SHOW(V) std::cout << #V " = " << V << std::endl
#define TRACKER_LOG(V)                                                                                               \
  std::cout << V;                                                                                                    \
  if (g_logfile)                                                                                                     \
  g_logfile << V

// ================================================================================================================
// STRUCTURES
// ================================================================================================================

struct relative_coord_type {
  int dx{0}, dy{0}; // Relative positions can be negative
};

struct subpix_coord_type {
  double x{0.0}, y{0.0};
};

struct imageDATA {
  float iso_speed;
  float shutter;
  float aperture;
  float focal_len;
  unsigned shot_order;
  std::string dateTime;
};

struct Grain {
  // variables initiales (generalement ce sont les donnees d'entree)
  int refcoord_xpix{0};   // coordonnée x de référence du point suivi (centre de masse)
  int refcoord_ypix{0};   // coordonnée y de référence du point suivi (centre de masse)
  double refrot{0.0};     // angle initial (de référence) du grain
  double radius_pix{0.0}; // Grain radius expressed in sub-pixels (input data)

  // Deplacements et rotation cumules par rapport aux variables initiales (par rapport à refxxx)
  double dx{0.0};
  double dy{0.0};
  double drot{0.0};

  // Backup des deplacements et rotations cumulés
  double dx_prev{0.0};
  double dy_prev{0.0};
  double drot_prev{0.0};

  // Kinematic data obtain on current DIC (num_image to num_image+iinc)
  double upix{0.0};    // x translation resulting from last DIC
  double vpix{0.0};    // y translation resulting from last DIC
  double rot_inc{0.0}; // Increment of rotation resulting from last DIC

  // Coefficients de correlation
  double NCC{0.0};              // Best NCC
  double NCC_rescue{0.0};       // Best NCC after rescue (0 if no rescue has been performed)
  double NCC_subpix{0.0};       // Best NCC after subpix (0 if no subpix refinement has been performed)
  double mean0{0.0}, C0C0{0.0}; // Precomputed data for the computation of NCCs
  bool masked{false};           // Pour les point dans une zone à ne pas correlé (utilisation avec une grille)

  std::vector<relative_coord_type> pattern;          // A list of relative positions that represent a pattern
  std::vector<relative_coord_type> pattern0_rotated; // Rotated pattern for the reference position
  std::vector<relative_coord_type> pattern1_rotated; // Rotated pattern for the tested position

  std::vector<int> neighbour; // List of neighbouring grains
  int num_neighbour{0};       // Number of neighbour for igrain (according to neighbour_dist_pix) @fixme supprimer

  Grain() = default;
  Grain(const Grain &c) = default;
  ~Grain() = default;

  void backup() {
    dx_prev = dx;
    dy_prev = dy;
    drot_prev = drot;
  }

  void reset() { NCC = NCC_rescue = NCC_subpix = 0.0; }
};

struct triangle {
  int i0{0}, i1{0}, i2{0};
};

struct increasingRadius_order {
  inline bool operator()(const Grain &G1, const Grain &G2) { return (G1.radius_pix < G2.radius_pix); }
};

struct decreasingRadius_order {
  inline bool operator()(const Grain &G1, const Grain &G2) { return (G1.radius_pix > G2.radius_pix); }
};

struct increasingHeight_order {
  inline bool operator()(const Grain &G1, const Grain &G2) { return (G1.refcoord_ypix < G2.refcoord_ypix); }
};

// Research zone (for pixel-precision tracking)
struct search_zone_type {
  // Be careful: they are positive values
  int left{1};         // Towards decreasing x
  int right{1};        // Towards increasing x
  int up{1};           // Towards decreasing y
  int down{1};         // Towards increasing y
  int num_rot{1};      // Number of tested angles (above and after, total = 2*num_rot+1)
  double inc_rot{1.0}; // Increment of rotation
};

struct data_stat {
  double ave;
  double adev;
  double sdev;
  double var;
  double skew;
  double curt;
};

// ================================================================================================================
// GLOBAL VARIABLES
// ================================================================================================================

// NOTE: Ces variables globales sont seulement DECLAREES ici (extern) afin que ce
// header puisse etre inclus par plusieurs unites de compilation (le CLI run.cpp et
// l'application tracker_gui.cpp), toutes deux liees a libtracker.a. Les DEFINITIONS
// (avec leurs valeurs par defaut) se trouvent dans tracker.cpp.

extern std::ofstream g_logfile;
extern int progress;

extern image_interpLinear IMAGE_INTERPOLATOR_LINEAR;
extern image_interpCubic IMAGE_INTERPOLATOR_CUBIC;
extern image_interpQuintic IMAGE_INTERPOLATOR_QUINTIC;
extern image_interpolator *IMAGE_INTERPOLATOR;

extern std::string procedure; // the default procedure is to track particles

extern int fake_undistor;

// Distortion parameters (grid method)
extern std::string grid_image_name;
extern int nx_grid_disto;
extern int ny_grid_disto;

// Distortion parameters (displacement method)
extern std::vector<int> image_numbers_corrDisto;
extern std::vector<double> imposed_displ_x_pix, imposed_displ_y_pix;
extern std::vector<std::vector<double>> dx_corrDisto, dy_corrDisto, NCC_subpix_corrDisto;
extern std::vector<double> disto_parameters;
extern std::vector<double> disto_parameters_perturb;
// 0 -> xc
// 1 -> yc
// 2 -> K1
// 3 -> K2
// 4 -> K3
// 5 -> P1
// 6 -> P2
// 7 -> P3

// These vectors are used for the correction of distortions
extern std::vector<int> List_Grains_i;
extern std::vector<int> List_Grains_j;
extern double equiProj_dist_min;
extern double equiProj_dist_max;

// subpixel positionning of rod centers
extern std::vector<double> xbound, ybound;
extern std::vector<double> center_parameters;
extern std::vector<double> center_parameters_perturb;
extern double subpix_center_dr0;
extern double subpix_center_xstep;
extern double subpix_center_ystep;
extern double subpix_center_rstep;
extern int subpix_center_threshold;

// this allows to define special pattern according to a specified radius
// A negative radius means "all the grains"
extern double targetRadiusPattern;

// Pattern quality
extern int halfPatternQual;

// Post-Process
extern int post_undistor;
extern int post_rescale;
extern int post_sync;
extern int post_rotation;
extern char device_data_name[256];
extern double radius_corners;
extern double radius_fixedPoints;
extern double scaling_distance;

extern int ShowSolidMotionError;

// Visualization
extern float colorMin;
extern float colorMax;
extern int visu_autoscale;
extern double visu_alpha;
extern std::string visu_output;
extern std::string visu_mode;
extern int nx_grid_visu;
extern int smoothDegreeF;
extern int visu_draw_in_ref;

// Minimization (Subpixel)

extern double subpix_tol;
extern double initial_direction_dx;
extern double initial_direction_dy;
extern double initial_direction_dr;

extern double subpix_displacement_max;
extern double overflow_stiffness;

// Core variables

extern int RawImages;

extern int DemosaicModel;
extern int rescaleGrayLevels;

extern std::vector<std::vector<std::vector<uint16_t>>> image; // image[0 or 1][dimx][dimy]

extern std::vector<imageDATA> imageData;                          // imageData[0 or 1]
extern char image_name[256];                                      // Image name in "printf" style.
extern char dic_name[256];                                        // DIC-file name in "printf" style.
extern int iref, ibeg, iend, iinc, iraz, idelta;                  // Control of file numbers
extern bool require_precomputations;                              // If true, remake the pre-computations
extern int dimx, dimy;                                            // Image sizes
extern std::string pattern_command;
extern std::string grain_positions_command;

extern int im_index_ref, im_index_current; // Index of image origin and target in table image
extern int wanted_num_threads;

extern int num_neighbour_max;
extern double neighbour_dist_pix; // max distance between points to find neighbours
extern int period_rebuild_neighbour_list;

// Flags
extern int verbose_level;
extern int rescue_level; // in range [0 2].
extern int subpixel;
extern int rotations;
extern int use_neighbour_list;

// Image generation for check
extern int make_images; // create images (1) or not (0)
extern double image_size; // in range [0 1] (DEPRECATED)
extern int image_div;   // will replace image_size (TODO)
extern int draw_angle;  //
extern int draw_disp;   //
extern int draw_rescued; //

extern std::vector<Grain> grain;
extern int num_grains;

extern std::vector<triangle> mesh;

extern std::vector<int> to_be_rescued;
extern int num_to_be_rescued;
extern double NCC_min; // The NCC-value above which a first rescue is mandatory

extern std::vector<int> to_be_super_rescued;
extern int num_to_be_super_rescued;
extern double NCC_min_super; // The NCC-value above which a super_rescue is mandatory

extern search_zone_type search_zone;              // For non sub-pixel tracking
extern search_zone_type search_zone_rescue;       // For points with bad mathing
extern search_zone_type search_zone_super_rescue; // For points with very bad mathing

// A table that hold the grain number being processed for each thread.
// It is required in subpixel tracking.
// Do not hesitate to increase the maximum number of threads if necessary
extern int igrain_of_thread[MAX_NUMBER_THREADS];

// ================================================================================================================
// Functions
// ================================================================================================================

// INTERACTION WITH USER  =========================================================================================
void header();
int init(int argc, char *argv[]);
void dialog();

// HELPERS ========================================================================================================
double get_time(); // mesure du temps (selon que OpenMP est utilisé ou pas)
inline int nearest(double x) { return (int)floor(x + 0.5); }
inline int nearest(float x) { return (int)floor(x + 0.5); }
inline double max(double a, double b) { return ((a > b) ? a : b); }
inline double max(float a, float b) { return ((a > b) ? a : b); }
void moment(std::vector<double> &data, data_stat &stat);

// CORE DIC PIT ===================================================================================================
void rotate_pixel_pattern(int igrain, int i, double c, double s, int *xpixel, int *ypixel);
void find_neighbours(int igrain);
void do_precomputations();
void follow_pattern_pixel(int igrain);
void follow_pattern_rescue_pixel(int igrain);
void follow_pattern_super_rescue_pixel(int igrain);
void follow_pattern_subpixel_xyR(int igrain);
double NCC_to_minimize_xyR(std::vector<double> &);
double Equiproj_Error_Computation(std::vector<double> &);
double disto_to_minimize_Equiproj(std::vector<double> &);
double disto_to_minimize_grid(std::vector<double> &);

// READ-WRITE FILES ===============================================================================================
int read_command_file(const char *name); // To read the command file
void write_data(const char *name);
void read_image(int i, const char *name, bool first_time);
void read_image(int i, int num, bool = false);
void read_raw_image(int i, const char *name, bool first_time);
int read_grains(const char *name, bool restart);
int read_gmsh(const char *name);
int save_grains(int num);
int save_grains(const char *name, int num, bool simpleVersion = false);
void create_image(int num);
void create_image_Netbpm(int num);

// PATTERNS =======================================================================================================
void make_grid(double xmin, double xmax, double ymin, double ymax, int nx, int ny, int aleaMax);
void make_grid_circular(double x_center, double y_center, double radius, double angle_inc);
void make_circ_pattern(int igrain, int radius_pix);
void make_ring_pattern(int igrain, int radius_IN_pix, int radius_OUT_pic);
void make_rect_pattern(int igrain, int half_width_pix, int half_height_pix);
void make_custom_pattern(const char *name);
void mask_rect(int xmin, int xmax, int ymin, int ymax);

// DISTO ==========================================================================================================
void distor(double *X, const double xd, const double yd, double &xu, double &yu);
void diff_distor(double *X, const double xu, const double yu, double &dxd_dxu, double &dyd_dyu);
double SolidMotionError_to_minimize(std::vector<double> &X);
void SolidMotionError(double &angle, double &xc, double &yc, double &xtrans, double &ytrans, double &mean,
                      double &stddev);
int undistor(double *X, const double xd, const double yd, double &xu, double &yu, int niterMax = 100,
             double tol = 1e-3, int nconv = 3);
void undistor_image(const char *name_from, const char *name_to);
void setSolidSolutionFromRef(double DX, double DY, double DROT);
void set_solidTranslation(double u, double v);

// PROCEDURES =====================================================================================================
void particle_tracking();
void particle_tracking_assisted_corrections();
void correction_distortion();
void correction_distortion_grid();
void find_subpixel_centers();
void pattern_quality();
void post_process();
void visu_process();
void gray_level_analysis();

// SYNTHETIC IMAGES ===============================================================================================
inline double findnoise2(double x, double y);
inline double interpolate(double a, double b, double x);
double noise(double x, double y);
void generate_synthetic_images(int w, int h, double zoom, int octaves, double p, double defx, double defy,
                               double transx_pix, double transy_pix, double rot_deg);

// ================================================================================================================