// g++-8 see.cpp -o see -I ~/Documents/devel/mbox/common  -framework OpenGL -framework GLUT -Wno-deprecated-declarations
// -DPNG_SCREENSHOTS 

#include "see.hpp"

void keyboardSpecial(int Key, int /*x*/, int /*y*/) {
  switch (Key) {
    case GLUT_KEY_UP:
      if (zslice > wslice) zslice -= wslice; 
      break;

    case GLUT_KEY_DOWN:
      if (zslice < data.dimz - wslice) zslice += wslice;
      break;
      
    case GLUT_KEY_LEFT:
      if (wslice > 2) wslice--; 
      break;

    case GLUT_KEY_RIGHT:
      if (wslice < 100) wslice++;
      break;
  };

  glutPostRedisplay();
}

void keyboard(unsigned char Key, int x, int y) {
  switch (Key) {
    case 27:  // ESCAPE
      //selectedParticle = -1;
      break;

    case ' ':
      textZone.reset();
      break;

    case 'b':
      show_background = 1 - show_background;
      break;

    case 'q':
      exit(0);
      break;

    case 's':
      show_slice = 1 - show_slice;
      break;


    case 'w': {
      vec3r d = center - eye;
      up.set(-d.x * d.y, d.x * d.x + d.z * d.z, -d.z * d.y);
      up.normalize();
    } break;


    case 'z': {
#ifdef PNG_H
      screenshot("oneshot.png");
#else
      screenshot("oneshot.tga");
#endif
    } break;
  
    case '=': {
      fitView();
      adjustClippingPlans();
    } break;
  };

  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    display();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
          mouse_mode = PAN;
        //else if (glutGetModifiers() == GLUT_ACTIVE_ALT)
        //  selection(x, y);
        else if (glutGetModifiers() == GLUT_ACTIVE_CTRL)
          mouse_mode = ZOOM;
        else
          mouse_mode = ROTATION;
        break;
      case GLUT_MIDDLE_BUTTON:
        mouse_mode = ZOOM;
        break;
    }
  }
}

void motion(int x, int y) {
  if (mouse_mode == NOTHING) return;

  double dx = (double)(x - mouse_start[0]) / (double)width;
  double dy = (double)(y - mouse_start[1]) / (double)height;
  double length;
  vec3r axis;

  switch (mouse_mode) {

    case ROTATION:
      axis = (cross(up, center - eye));
      axis.normalize();
      eye = geoTool::rotatePoint(eye, center, up, -dx * M_PI);
      eye = geoTool::rotatePoint(eye, center, axis, dy * M_PI);
      up = (geoTool::rotatePoint((center + up), center, axis, dy * M_PI) - center);
      up.normalize();
      break;

    case ZOOM:
      eye = center + (eye - center) * (dy + 1.0);
      break;

    case PAN:
      length = (eye - center).length() * tan(view_angle * M_PI / 360.0) * 2.0;
      axis = cross(up, center - eye);
      axis.normalize();
      center = center + axis * dx * length * 0.8;
      center = center + up * dy * length;
      break;

    default:
      break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  display();
}

void display() {
  glutTools::clearBackground(show_background);
  adjustClippingPlans();
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  if (show_slice) drawScan();
  if (show_shells) drawShells();
  if (show_zone) drawZone();
  drawSizes();

  textZone.draw();
  glFlush();
  glutSwapBuffers();
}

void computePerspective() {
  double zf = (eye - center).normalize();

  vec3r mx(data.dimx, data.dimz, data.dimy);
  max_length = (GLfloat)(2.0 * norm(mx));

  znear = zf - 0.5 * max_length;
  double close_dst = 0.1 * zf;
  if (znear < close_dst) znear = close_dst;
  zfar = zf + 0.5 * max_length;
}

void adjustClippingPlans() {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  wh_ratio = (float)width / (float)height;
  computePerspective();
  gluPerspective(view_angle, wh_ratio, znear, zfar);
  glMatrixMode(GL_MODELVIEW);
}

void fitView() {
  vec3r dir = (eye - center);
  vec3r diag(data.dimx, -(double)data.dimz, data.dimy);
  dir.normalize();
  center = 0.5 * diag;
  GLfloat d = 0.5 * diag.length() / (atan(view_angle * M_PI / 360.0));
  eye = center + d * dir;
}

void reshape(int w, int h) {
  width = w;
  height = h;
  glViewport(0, 0, width, height);

  adjustClippingPlans();
  glutPostRedisplay();
}


void drawSizes() {
  glLineWidth(1.0f);
  glColor4ub(255, 0, 0, 255);  // RED
  
  vec3r corner;
  glBegin(GL_LINE_LOOP);
  glVertex3i(0, 0, 0);
  glVertex3i(0, 0         , data.dimy);
  glVertex3i(0, -data.dimz, data.dimy);
  glVertex3i(0, -data.dimz, 0        );
  glEnd();
  glBegin(GL_LINE_LOOP);
  glVertex3i(data.dimx, 0         , 0        );
  glVertex3i(data.dimx, 0         , data.dimy);
  glVertex3i(data.dimx, -data.dimz, data.dimy);
  glVertex3i(data.dimx, -data.dimz, 0        );
  glEnd();
  
  glBegin(GL_LINES);
  glVertex3i(0, 0, 0);
  glVertex3i(data.dimx, 0         , 0        );
  
  glVertex3i(0, 0         , data.dimy);
  glVertex3i(data.dimx, 0         , data.dimy);
  
  glVertex3i(0, -data.dimz, data.dimy);
  glVertex3i(data.dimx, -data.dimz, data.dimy);
  
  glVertex3i(0, -data.dimz, 0        );
  glVertex3i(data.dimx, -data.dimz, 0        );
  glEnd();
  
}

void drawShells() {
  
  glCullFace(GL_FRONT_AND_BACK);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glColor3ub(0, 250, 0);
  glPointSize(2.0);
  glBegin(GL_POINTS);
  
  int nAngle = data.nbAngle;
  int nH = data.nbH;
  
  for (size_t ish = 0 ; ish < data.shells.size() ; ish++) {

    double dAngle = 2.0 * M_PI / (double)nAngle;
    double dh = data.shells[ish].h / (double)nH;
  
    for (double angle = 0.0 ; angle < 2.0 * M_PI ; angle += dAngle) {
      double yin = data.shells[ish].r_in * cos(angle);
      double zin = data.shells[ish].r_in * sin(angle);
      double yout = data.shells[ish].r_out * cos(angle);
      double zout = data.shells[ish].r_out * sin(angle);
      for (double x = -0.5 * data.shells[ish].h ; x <= 0.5 * data.shells[ish].h; x += dh) {
        double vx = x;
        double vy = yout;
        double vz = zout;
        rotyz(data.shells[ish].roty, data.shells[ish].rotz, vx, vy, vz);
        vx += data.shells[ish].x;
        vy += data.shells[ish].y;
        vz += data.shells[ish].z;
        glVertex3f (vx, -vz, vy);
        
        vx = x;
        vy = yin;
        vz = zin;
        rotyz(data.shells[ish].roty, data.shells[ish].rotz, vx, vy, vz);
        vx += data.shells[ish].x;
        vy += data.shells[ish].y;
        vz += data.shells[ish].z;
        glVertex3f (vx, -vz, vy);
      }
    }
  }

  glEnd();
}

void drawScan() {
  if (mouse_mode != NOTHING) {
    drawSizes();
    return;
  }
  
  glCullFace(GL_FRONT_AND_BACK);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glColor3ub(178, 34, 34);
  glPointSize(2.0);
  glBegin(GL_POINTS);
  size_t id = 0;
  size_t nxny = data.dimx * data.dimy;
  size_t nx = data.dimx;
  size_t zmin = zslice;
  if (zslice > wslice) zmin = zslice - wslice;
  size_t zmax = zslice;
  if (zslice < data.dimz - wslice) zmax = zslice + wslice;
  for (size_t z = zmin ; z < zmax ; z += 1) {
    for (size_t y = 0 ; y < data.dimy ; y++) {
      for (size_t x = 0 ; x < data.dimx ; x++) {
        id = z * nxny + y * nx + x;
        if (data.image[id] == 1) {
          glVertex3i (x, -z, y);
        }
        id++;
      }
    }
  }
  glEnd();
}


void drawZone() {
  glCullFace(GL_FRONT_AND_BACK);
  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
  glColor3ub(0, 0, 250);
  glPointSize(2.0);
  
  glBegin(GL_LINES);
  
  int nAngle = 18;
  double dAngle = 2.0 * M_PI / (double)nAngle;
  
  for (double angle = 0.0 ; angle < 2.0 * M_PI ; angle += dAngle) {
    double ry0 = data.zone.R0 * cos(angle);
    double rx0 = data.zone.R0 * sin(angle);
    double ry1 = data.zone.R1 * cos(angle);
    double rx1 = data.zone.R1 * sin(angle);
        
    glVertex3f (data.zone.x0 + rx0, -data.zone.z0, data.zone.y0 +ry0);
    glVertex3f (data.zone.x1 + rx1, -data.zone.z1, data.zone.y1 +ry1);

  }
  glEnd();
}

// This function is for making a screenshot (or a series of screenshots)
// The output format is PNG if linked with libpng, else it is TGA.
int screenshot(const char* filename) {
#ifdef PNG_H
  int i;

  int screenStats[4];
  glGetIntegerv(GL_VIEWPORT, screenStats);
  int x = screenStats[0];
  int y = screenStats[1];
  int width = screenStats[2];
  int height = screenStats[3];

  FILE* fp;
  png_structp png_ptr;
  png_infop info_ptr;
  png_bytep* rowp;
  GLubyte* glimage = NULL;

  fp = fopen(filename, "wb");
  if (!fp) return 1;

  glimage = (GLubyte*)malloc(width * height * sizeof(GLubyte) * 3);

  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glReadPixels(x, y, width, height, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)glimage);

  rowp = (png_bytep*)malloc(sizeof(png_bytep*) * height);
  if (!rowp) return 1;

  for (i = 0; i < height; i++) {
    rowp[i] = (png_bytep)&glimage[3 * ((height - i - 1) * width)];
  }

  png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) return 1;

  info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) return 1;

  png_init_io(png_ptr, fp);
  png_set_IHDR(png_ptr, info_ptr, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
               PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
  png_write_info(png_ptr, info_ptr);
  png_write_image(png_ptr, rowp);
  png_write_end(png_ptr, info_ptr);
  png_destroy_write_struct(&png_ptr, &info_ptr);

  free(rowp);
  fflush(stdout);
  fclose(fp);
  return 0;

#else

  // http://forum.devmaster.net/t/rendering-a-single-frame-to-a-file-with-opengl/12469/2

  // we will store the image data here
  unsigned char* pixels;
  // the thingy we use to write files
  FILE* shot;
  // we get the width/height of the screen into this array
  int screenStats[4];

  // get the width/height of the window
  glGetIntegerv(GL_VIEWPORT, screenStats);

  // generate an array large enough to hold the pixel data
  // (width*height*bytesPerPixel)
  pixels = new unsigned char[screenStats[2] * screenStats[3] * 3];
  // read in the pixel data, TGA's pixels are BGR aligned
  glReadPixels(0, 0, screenStats[2], screenStats[3], GL_BGR, GL_UNSIGNED_BYTE, pixels);

  // open the file for writing. If unsucessful, return 1
  if ((shot = fopen(filename, "wb")) == NULL) return 1;

  // this is the tga header it must be in the beginning of
  // every (uncompressed) .tga
  unsigned char TGAheader[12] = {0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  // the header that is used to get the dimensions of the .tga
  // header[1]*256+header[0] - width
  // header[3]*256+header[2] - height
  // header[4] - bits per pixel
  // header[5] - ?
  unsigned char header[6] = {((unsigned char)(screenStats[2] % 256)),
                             ((unsigned char)(screenStats[2] / 256)),
                             ((unsigned char)(screenStats[3] % 256)),
                             ((unsigned char)(screenStats[3] / 256)),
                             24,
                             0};

  // write out the TGA header
  fwrite(TGAheader, sizeof(unsigned char), 12, shot);
  // write out the header
  fwrite(header, sizeof(unsigned char), 6, shot);
  // write the pixels
  fwrite(pixels, sizeof(unsigned char), screenStats[2] * screenStats[3] * 3, shot);

  // close the file
  fclose(shot);
  // free the memory
  delete[] pixels;

  // return success
  return 0;
#endif
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {
  
  if (argc == 2) {
    read_data(argv[1], data);
  }
  else {
    std::cout << "Usage: "<< argv[0] << " inputData.txt\n\n";
    exit(0);
  }
  
  // READ the raw image
  char fname[256];
  sprintf(fname, "%s.raw", data.filename.c_str());
  std::cout << "Loading rawFile " << fname << std::endl;
  data.image = readCompressedRaw(fname, data.dimx, data.dimy, data.dimz);
  zslice = data.dimz / 2;
  
  // READ the shells
  sprintf(fname, "%s.txt", data.filename.c_str());
  std::cout << "Loading shellFile " << fname << std::endl; /// <<<<<
  readShells(fname, data.shells);  
  
  // display info
  textZone.addLine("filename: %s (.raw, .txt)", data.filename.c_str());
  textZone.addLine("size: %zu x %zu x %zu", data.dimx, data.dimy, data.dimz);

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("Rockable VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(keyboardSpecial);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Init the visualizer
  center.set(0.0, 0.0, 0.0);  // where we look at
  eye.set(0.0, 0.0, 1.0);     // from where we look
  up.set(0.0, 1.0, 0.0);      // direction (normalized)

  mouse_mode = NOTHING;
  view_angle = 45.0;
  znear = 0.01;
  zfar = 10.0;

  glDisable(GL_CULL_FACE);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);
  glEnable(GL_COLOR_MATERIAL);

  // Create light components
  GLfloat ambientLight[] = {0.2f, 0.2f, 0.2f, 1.0f};
  GLfloat diffuseLight[] = {0.8f, 0.8f, 0.8, 1.0f};
  GLfloat specularLight[] = {0.5f, 0.5f, 0.5f, 1.0f};
  GLfloat positionLight0[] = {1000000.0f, 1000000.0f, 1000000.0f, 1.0f};
  GLfloat positionLight1[] = {-1000000.0f, -1000000.0f, -1000000.0f, 1.0f};

  // Assign created components to GL_LIGHT0
  glLightfv(GL_LIGHT0, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT0, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT0, GL_POSITION, positionLight0);

  // Assign created components to GL_LIGHT1
  glLightfv(GL_LIGHT1, GL_AMBIENT, ambientLight);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuseLight);
  glLightfv(GL_LIGHT1, GL_SPECULAR, specularLight);
  glLightfv(GL_LIGHT1, GL_POSITION, positionLight1);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_POINT_SMOOTH);
  //glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LEQUAL);

  // ==== Enter GLUT event processing cycle
  adjustClippingPlans();
  fitView();
  glutMainLoop();
  return 0;
}
