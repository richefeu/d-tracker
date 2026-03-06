#ifndef SEE_DIC_HPP
#define SEE_DIC_HPP

#include <GL/freeglut.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>

#include "dicProcess.hpp"
#include "ColorTable.hpp"

DICResult dic_out;
int dicNum = 1;

AABB viewBox;

int main_window;

// flags
int show_background = 1;
int show_particles = 1;
int show_velocities = 0;
int show_cell = 1;
int show_forces = 0;
int showOrientations = 0;

int width = 1032;
int height = 776;
float wh_ratio = (float)width / (float)height;

float colNCCMin = 0.7;
int colorMode = 1;
int iSelected = -1;

// Miscellaneous global variables
enum MouseMode { NOTHING, ROTATION, ZOOM, PAN } mouse_mode = NOTHING;
int display_mode = 0;  // sample or slice rotation
int mouse_start[2];

// Drawing functions
void drawParticlesQuickly();
void drawParticles();
void print(void *font, int x, int y, const char *fmt, ...);
void drawColorMap();


// Callback functions
void keyboard(unsigned char Key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void display();
void reshape(int x, int y);

// Helper functions
void printHelp();
void fit_view();
void findSelected(int x, int y);
bool fileExists(const char* fileName);
void try_to_readDIC(int num);

#endif /* end of include guard: SEE_CONF_HPP */
