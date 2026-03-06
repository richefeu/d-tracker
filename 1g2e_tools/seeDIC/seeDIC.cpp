#include "seeDIC.hpp"
#include "dicProcess.cpp"

void printHelp() {
  using namespace std;
  cout << endl;
  cout << "+         load next DIC file" << endl;
  cout << "-         load previous DIC file" << endl;
  cout << "=         fit the view" << endl;
  cout << "0/1/2/3   colors of grains (all grey, NCC, NCC_rescue, NCC_subpix)" << endl;
  cout << "space     select the grain under the mouse" << endl;
  cout << "i         print data of selected grain" << endl;
  cout << "n/N       tune the lower value of NCC" << endl;
  cout << "o         show/hide rotations" << endl;
  cout << "g         go to a console-specified file number" << endl;
  cout << "q         quit" << endl;
  // cout << "" << endl;
  cout << endl;
}

void keyboard(unsigned char Key, int x, int y) {
  switch (Key) {

    case 'q': {
      exit(0);
    } break;

    case 'g': {
      cout << "Go to file number: ";
      int conNumTry;
      cin >> conNumTry;
      try_to_readDIC(conNumTry);
    } break;

    case 'o': {
      showOrientations = 1 - showOrientations;
    } break;

    case 'n': {
      colNCCMin -= 0.025;
      if (colNCCMin < 0.0) colNCCMin = 0.0;
      std::cout << "colNCCMin = " << colNCCMin << '\n';
    } break;
    case 'N': {
      colNCCMin += 0.025;
      if (colNCCMin > 0.975) colNCCMin = 0.975;
      std::cout << "colNCCMin = " << colNCCMin << '\n';
    } break;

    case '0': {
      colorMode = 0;
      std::cout << "no colors\n";
    } break;
    case '1': {
      colorMode = 1;
      std::cout << "colored by NCC\n";
    } break;
    case '2': {
      colorMode = 2;
      std::cout << "colored by NCC_rescue\n";
    } break;
    case '3': {
      colorMode = 3;
      std::cout << "colored by NCC_subpix\n";
    } break;

    case ' ': {
      findSelected(x, y);
    } break;

    case 'i':
    case 'I': {
      if (iSelected > 0 && (size_t)iSelected < dic_out.grains.size()) {
        std::cout << "===================================================================\n";
        std::cout << "iSelected     = " << iSelected << ", line " << 2 + iSelected << " in the file dic_out_" << dicNum
                  << ".txt\n";
        std::cout << "refcoord_xpix = " << dic_out.grains[iSelected].refcoord_xpix << '\n';
        std::cout << "refcoord_ypix = " << dic_out.grains[iSelected].refcoord_ypix << '\n';
        std::cout << "refrot        = " << dic_out.grains[iSelected].refrot << '\n';
        std::cout << "radius_pix    = " << dic_out.grains[iSelected].radius_pix << '\n';
        std::cout << "dx            = " << dic_out.grains[iSelected].dx << '\n';
        std::cout << "dy            = " << dic_out.grains[iSelected].dy << '\n';
        std::cout << "drot          = " << dic_out.grains[iSelected].drot << '\n';
        std::cout << "upix          = " << dic_out.grains[iSelected].upix << '\n';
        std::cout << "vpix          = " << dic_out.grains[iSelected].vpix << '\n';
        std::cout << "rot_inc       = " << dic_out.grains[iSelected].rot_inc << '\n';
        std::cout << "NCC           = " << dic_out.grains[iSelected].NCC << '\n';
        std::cout << "NCC_rescue    = " << dic_out.grains[iSelected].NCC_rescue << '\n';
        std::cout << "NCC_subpix    = " << dic_out.grains[iSelected].NCC_subpix << '\n';
        std::cout << "===================================================================\n";
      }

    } break;

    case '-': {
      if (dicNum > 0) try_to_readDIC(dicNum - 1);
    } break;

    case '+': {
      try_to_readDIC(dicNum + 1);
    } break;

    case '=': {
      fit_view();
    } break;
  };

  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    // display();
    glutPostRedisplay();
  } else if (state == GLUT_DOWN) {
    mouse_start[0] = x;
    mouse_start[1] = y;
    switch (button) {
      case GLUT_LEFT_BUTTON:
        if (glutGetModifiers() == GLUT_ACTIVE_SHIFT)
          mouse_mode = PAN;
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

  switch (mouse_mode) {

    case ZOOM: {
      double ddy = (viewBox.ymax - viewBox.ymin) * dy;
      double ddx = (viewBox.xmax - viewBox.xmin) * dy;
      viewBox.xmin -= ddx;
      viewBox.xmax += ddx;
      viewBox.ymin -= ddy;
      viewBox.ymax += ddy;
    } break;

    case PAN: {
      double ddx = (viewBox.xmax - viewBox.xmin) * dx;
      double ddy = (viewBox.ymax - viewBox.ymin) * dy;
      viewBox.xmin -= ddx;
      viewBox.xmax -= ddx;
      viewBox.ymin -= ddy;
      viewBox.ymax -= ddy;
    } break;

    default:
      break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  reshape(width, height);
  // display();
  glutPostRedisplay();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  if (mouse_mode == NOTHING) {
    drawParticles();
  } else {
    drawParticlesQuickly();
  }

  if (colorMode > 0) drawColorMap();

  glFlush();
  glutSwapBuffers();
}

void findSelected(int x, int y) {

  // copy-past of what we've got in reshape --------
  GLfloat aspect = (GLfloat)width / (GLfloat)height;
  double left = viewBox.xmin;
  double right = viewBox.xmax;
  double bottom = viewBox.ymin;
  double top = viewBox.ymax;
  double worldW = right - left;
  double worldH = top - bottom;
  double dW = 0.1 * worldW;
  double dH = 0.1 * worldH;
  left -= dW;
  right += dW;
  top += dH;
  bottom -= dH;
  worldW = right - left;
  worldH = top - bottom;
  if (worldW >= worldH) {
    worldH = worldW / aspect;
    top = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * aspect;
    right = 0.5 * (left + right + worldW);
    left = right - worldW;
  }
  // -----------------------------------------------

  double xworld = left + (x / (double)width) * (right - left);
  double yworld = bottom + (y / (double)height) * (top - bottom);

  iSelected = -1;
  for (size_t i = 0; i < dic_out.grains.size(); ++i) {
    double xc = dic_out.grains[i].refcoord_xpix + dic_out.grains[i].dx;
    double yc = dic_out.grains[i].refcoord_ypix + dic_out.grains[i].dy;
    double R = dic_out.grains[i].radius_pix;

    double dx = xc - xworld;
    double dy = yc - yworld;
    double dst2 = dx * dx + dy * dy;
    if (dst2 < R * R) {
      iSelected = i;
      break;
    }
  }
}

void fit_view() {
  viewBox.xmin = viewBox.xmax = dic_out.grains[0].refcoord_xpix + dic_out.grains[0].dx;
  viewBox.ymin = viewBox.ymax = dic_out.grains[0].refcoord_ypix + dic_out.grains[0].dy;
  double R0 = dic_out.grains[0].radius_pix;
  viewBox.xmin -= R0;
  viewBox.xmax += R0;
  viewBox.ymin -= R0;
  viewBox.ymax += R0;
  for (size_t i = 1; i < dic_out.grains.size(); ++i) {
    double xc = dic_out.grains[i].refcoord_xpix + dic_out.grains[i].dx;
    double yc = dic_out.grains[i].refcoord_ypix + dic_out.grains[i].dy;
    double R = dic_out.grains[i].radius_pix;
    double xmin = xc - R;
    double xmax = xc + R;
    double ymin = yc - R;
    double ymax = yc + R;
    if (xmin < viewBox.xmin) viewBox.xmin = xmin;
    if (xmax > viewBox.xmax) viewBox.xmax = xmax;
    if (ymin < viewBox.ymin) viewBox.ymin = ymin;
    if (ymax > viewBox.ymax) viewBox.ymax = ymax;
  }

  reshape(width, height);
}

void reshape(int w, int h) {
  width = w;
  height = h;
  GLfloat aspect = (GLfloat)width / (GLfloat)height;
  double left = viewBox.xmin;
  double right = viewBox.xmax;
  double bottom = viewBox.ymin;
  double top = viewBox.ymax;
  double worldW = right - left;
  double worldH = top - bottom;
  double dW = 0.1 * worldW;
  double dH = 0.1 * worldH;
  left -= dW;
  right += dW;
  top += dH;
  bottom -= dH;
  worldW = right - left;
  worldH = top - bottom;
  if (worldW >= worldH) {
    worldH = worldW / aspect;
    top = 0.5 * (bottom + top + worldH);
    bottom = top - worldH;
  } else {
    worldW = worldH * aspect;
    right = 0.5 * (left + right + worldW);
    left = right - worldW;
  }

  glViewport(0, 0, width, height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(left, right, top, bottom);  // inverse top and bottom because y-axis is toward bottom
}

void print(void* font, int x, int y, const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);
  char buffer[128];
  vsnprintf(buffer, 127, fmt, args);
  va_end(args);

  glRasterPos2i(x, y);
  for (size_t i = 0; buffer[i]; ++i) glutBitmapCharacter(font, buffer[i]);
}

void drawColorMap() {
  // switch2D
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, width, 0, height, -1.0f, 1.0f);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glPushAttrib(GL_DEPTH_TEST);
  glDisable(GL_DEPTH_TEST);
  glPushAttrib(GL_LIGHTING);
  glDisable(GL_LIGHTING);

  glColor3i(0, 0, 0);
  switch (colorMode) {
    case 1:
      print(GLUT_BITMAP_8_BY_13, 10, 10, "NCC");
      break;

    case 2:
      print(GLUT_BITMAP_8_BY_13, 10, 10, "NCC_rescue");
      break;

    case 3:
      print(GLUT_BITMAP_8_BY_13, 10, 10, "NCC_subpix");
      break;

    default:
      print(GLUT_BITMAP_8_BY_13, 10, 10, "???");
      break;
  }

  ColorTable CT;
  CT.setSize(16);
  CT.setSwap(true);
  CT.setMinMax(colNCCMin, 1.0);
  CT.Rebuild();
  colorRGBA col;

  float res = 0.01;
  float h = 0.8 * height;
  for (float val = 0.0; val <= 1.0; val += res) {
    CT.getRGB(val, &col);
    glColor3ub(col.r, col.g, col.b);
    glBegin(GL_POLYGON);
    glVertex2f(5, 30 + val * h);
    glVertex2f(20, 30 + val * h);
    glVertex2f(20, 30 + (val + res) * h);
    glVertex2f(5, 30 + (val + res) * h);
    glEnd();
  }

  glColor3ub(0, 0, 0);
  print(GLUT_BITMAP_8_BY_13, 25, 30, "0.0");
  print(GLUT_BITMAP_8_BY_13, 25, 30 + h, "1.0");
  print(GLUT_BITMAP_8_BY_13, 25, 30 + h * colNCCMin, "%0.2f", colNCCMin);

  // switch3D back
  glPopAttrib();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}

void drawParticlesQuickly() {

  glLineWidth(1.0f);
  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
  for (size_t i = 0; i < dic_out.grains.size(); ++i) {
    double xc = dic_out.grains[i].refcoord_xpix + dic_out.grains[i].dx;
    double yc = dic_out.grains[i].refcoord_ypix + dic_out.grains[i].dy;
    double R = dic_out.grains[i].radius_pix;
    glBegin(GL_LINE_LOOP);
    for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.1 * M_PI) {
      glVertex2f(xc + R * cos(angle), yc + R * sin(angle));
    }
    glEnd();
  }
}

void drawParticles() {
  ColorTable CT;
  CT.setSize(16);
  CT.setSwap(true);
  CT.setMinMax(colNCCMin, 1.0);
  CT.Rebuild();
  colorRGBA col;

  for (size_t i = 0; i < dic_out.grains.size(); ++i) {
    double xc = dic_out.grains[i].refcoord_xpix + dic_out.grains[i].dx;
    double yc = dic_out.grains[i].refcoord_ypix + dic_out.grains[i].dy;
    double R = dic_out.grains[i].radius_pix;

    if (R == 1.0) {  // Corners
      glLineWidth(1.0f);
      glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
      glBegin(GL_LINES);
      glVertex2f(xc - 40, yc - 40);
      glVertex2f(xc + 40, yc + 40);
      glVertex2f(xc - 40, yc + 40);
      glVertex2f(xc + 40, yc - 40);
      glEnd();
    } else if (R == 2.1) {  // scaling points
      glLineWidth(1.0f);
      glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
      glBegin(GL_LINES);
      glVertex2f(xc, yc - 40);
      glVertex2f(xc, yc + 40);
      glVertex2f(xc - 40, yc);
      glVertex2f(xc + 40, yc);
      glEnd();
    } else if (R == 2.2) {  // fixe points
      glPointSize(10.0f);
      glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
      glBegin(GL_POINTS);
      glVertex2f(xc, yc);
      glEnd();
    } else {
      if ((int)i == iSelected) {
        glLineWidth(4.0f);
      } else {
        glLineWidth(1.0f);
      }

      switch (colorMode) {
        case 1: {
          CT.getRGB(dic_out.grains[i].NCC, &col);
          glColor3ub(col.r, col.g, col.b);
        } break;
        case 2: {
          CT.getRGB(dic_out.grains[i].NCC_rescue, &col);
          glColor3ub(col.r, col.g, col.b);
        } break;
        case 3: {
          CT.getRGB(dic_out.grains[i].NCC_subpix, &col);
          glColor3ub(col.r, col.g, col.b);
        } break;
        default: {
          glColor3f(0.9f, 0.9f, 0.9f);
        }
      }

      glBegin(GL_POLYGON);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(xc + R * cos(angle), yc + R * sin(angle));
      }
      glEnd();

      glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
      glBegin(GL_LINE_LOOP);
      for (double angle = 0.0; angle < 2.0 * M_PI; angle += 0.05 * M_PI) {
        glVertex2f(xc + R * cos(angle), yc + R * sin(angle));
      }
      glEnd();

      if (showOrientations) {
        double drot = dic_out.grains[i].refrot + dic_out.grains[i].drot;
        glBegin(GL_LINES);
        glVertex2f(xc, yc);
        glVertex2f(xc + R * cos(drot), yc + R * sin(drot));
        glEnd();
      }
    }
  }
}

/// Robust and portable function to test if a file exists
bool fileExists(const char* fileName) {
  std::fstream fin;
  fin.open(fileName, std::ios::in);
  if (fin.is_open()) {
    fin.close();
    return true;
  }
  fin.close();
  return false;
}

void try_to_readDIC(int num) {
  char file_name[256];
  sprintf(file_name, "dic_out_%d.txt", num);
  if (fileExists(file_name)) {
    std::cout << "Read " << file_name << std::endl;
    dicNum = num;
    dic_out.read(".", num);
  } else
    std::cout << file_name << " does not exist" << std::endl;
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {

  if (argc == 1) {
    try_to_readDIC(1);
  } else if (argc == 2) {
    try_to_readDIC(atoi(argv[1]));
  } else if (argc == 3) {
    try_to_readDIC(atoi(argv[1]));
    // try_to_readCommands(atoi(argv[2]));
  }

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("DIC VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  // buildMenu();
  // glutAttachMenu(GLUT_RIGHT_BUTTON);

  mouse_mode = NOTHING;

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
