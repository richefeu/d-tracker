#ifndef SEE_HPP_E29BD15E
#define SEE_HPP_E29BD15E

# ifdef __APPLE__
#  include <GLUT/glut.h>
# else
#  include <GL/glut.h>
# endif

# ifdef PNG_SCREENSHOTS
#  include <png.h>
# endif

#include <vector>

#include "glutTools.hpp"
#include "geoTool.hpp"

#include "rawFiles.hpp"
#include "shellFiles.hpp"

Data data;

int main_window;

// flags with default values
int show_background    = 0;
int show_slice         = 1;
int show_shells        = 1;
int show_zone          = 1;
int show_sliced_shells = 1;
int show_keybinds      = 0;

GLfloat alpha_particles = 1.0f;
GLfloat alpha_fixparticles = 0.1f;

int width  = 800;
int height = 800;
float wh_ratio = (float)width / (float)height;

glTextZone textZone(3, &width, &height);

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int mouse_start[2];
float view_angle;
float znear;
float zfar;
GLfloat Rot_Matrix[16];
GLfloat max_length;

vec3r eye;
vec3r center;
vec3r up;

size_t zslice = 0; // z position of slice
size_t wslice = 10; // width of slice

void drawSizes();
void drawScan();
void drawShells();
void drawZone();

// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion (int x, int y);
void display();
void reshape(int x, int y);

// Helper functions
void computePerspective();
void adjustClippingPlans();
void fitView();
int screenshot(const char *filename);

#endif /* end of include guard: SEE_HPP_E29BD15E */
