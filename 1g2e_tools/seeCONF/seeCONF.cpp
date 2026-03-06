#include "seeCONF.hpp"
#include "confProcess.cpp"

void printHelp() {
  using namespace std;
  cout << endl;
  cout << "+         load next configuration file" << endl;
  cout << "-         load previous configuration file" << endl;
  cout << "=         fit the view" << endl;
  cout << "q         quit" << endl;
  // cout << "" << endl;
  cout << endl;
}

void printInfo() {
  using namespace std;

  cout << "Reference Conf = " << refConfNum << "\n";
  cout << "Current Conf = " << confNum << "\n";

  cout << "Box:\n";
  cout << "top left = (" << Conf.Box.x[0] << ", " << Conf.Box.y[0] << ")\n";
  cout << "top right = (" << Conf.Box.x[1] << ", " << Conf.Box.y[1] << ")\n";
  cout << "bottom right = (" << Conf.Box.x[2] << ", " << Conf.Box.y[2] << ")\n";
  cout << "bottom left = (" << Conf.Box.x[3] << ", " << Conf.Box.y[3] << ")\n";
}

void keyboard(unsigned char Key, int x, int y) {
  switch (Key) {

    case '0': {
      color_option = 0;
    } break;

    case '1': {
      colorTable.setMinMax(0.5, 1.0);
      colorTable.setTableID(2);
      colorTable.Rebuild();
      color_option = 1;
    } break;

    case '2': {
      colorTable.setMinMax(0.5, 1.0);
      colorTable.setTableID(2);
      colorTable.Rebuild();
      color_option = 2;
    } break;

    case 'i': {
      printInfo();
    } break;

    case 'S': {
      vScale *= 1.05;
    } break;

    case 's': {
      vScale *= 0.95;
      if (vScale < 0.0) vScale = 1.0;
    } break;

    case 'f': {
      show_fluctuations = 1 - show_fluctuations;
      if (show_fluctuations == 1 && show_displacements == 1) show_displacements = 0;
    } break;
    
    case 'v': {
      show_displacements = 1 - show_displacements;
      if (show_displacements == 1 && show_fluctuations == 1) show_fluctuations = 0;
    } break;

    case 'q': {
      exit(0);
    } break;

    case 'g': {
      cout << "Go to file number ";
      int conNumTry;
      cin >> conNumTry;
      try_to_readConf(conNumTry, Conf, confNum);
    } break;

    case 'r': {
      ref_fixed = 1 - ref_fixed;
    } break;

    case '-': {
      if (ref_fixed == 0 && refConfNum > 1) {
        std::cout << "Reference Configuration: ";
        try_to_readConf(refConfNum - 1, RefConf, refConfNum);
      }
      std::cout << "Current Configuration: ";
      if (confNum > 0) try_to_readConf(confNum - 1, Conf, confNum);
    } break;

    case '+': {
      if (ref_fixed == 0) {
        std::cout << "Reference Configuration: ";
        try_to_readConf(refConfNum + 1, RefConf, refConfNum);
      }
      std::cout << "Current Configuration: ";
      try_to_readConf(confNum + 1, Conf, confNum);
    } break;

    case '=': {
      fit_view();
      reshape(width, height);
    } break;
  };

  glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {

  if (state == GLUT_UP) {
    mouse_mode = NOTHING;
    //display();
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
      double ddy = (worldBox.ymax - worldBox.ymin) * dy;
      double ddx = (worldBox.xmax - worldBox.xmin) * dy;
      worldBox.xmin -= ddx;
      worldBox.xmax += ddx;
      worldBox.ymin -= ddy;
      worldBox.ymax += ddy;
    } break;

    case PAN: {
      double ddx = (worldBox.xmax - worldBox.xmin) * dx;
      double ddy = (worldBox.ymax - worldBox.ymin) * dy;
      worldBox.xmin -= ddx;
      worldBox.xmax -= ddx;
      worldBox.ymin += ddy;
      worldBox.ymax += ddy;
    } break;

    default:
      break;
  }
  mouse_start[0] = x;
  mouse_start[1] = y;

  reshape(width, height);
  //display();
  glutPostRedisplay();
}

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);

  drawBox();
  drawParticles();
  if (show_displacements == 1) drawDisplacements();
  if (show_fluctuations == 1) drawFluctuations();

  glFlush();
  glutSwapBuffers();
}

void fit_view() { worldBox = Conf.Box.getAABB(); }

void reshape(int w, int h) {
  width = w;
  height = h;
  GLfloat aspect = (GLfloat)width / (GLfloat)height;
  double left = worldBox.xmin;
  double right = worldBox.xmax;
  double bottom = worldBox.ymin;
  double top = worldBox.ymax;
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
  gluOrtho2D(left, right, bottom, top);

  glutPostRedisplay();
}

void drawBox() {
  glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
  glLineWidth(2.0f);

  glBegin(GL_LINE_LOOP);
  glVertex2f(Conf.Box.x[0], Conf.Box.y[0]);
  glVertex2f(Conf.Box.x[1], Conf.Box.y[1]);
  glVertex2f(Conf.Box.x[2], Conf.Box.y[2]);
  glVertex2f(Conf.Box.x[3], Conf.Box.y[3]);
  glEnd();

  double Rcorner = 0.02;  // 2 cm is what we get on the device 1g2e
  for (int cc = 0; cc < 4; cc++) {
    glBegin(GL_LINE_STRIP);
    double angle0 = 0.5 * M_PI * (1 + cc) - 0.25 * M_PI;
    for (double angle = angle0; angle < angle0 + M_PI; angle += 0.05 * M_PI) {
      glVertex2f(Conf.Box.x[cc] + Rcorner * cos(angle), Conf.Box.y[cc] + Rcorner * sin(angle));
    }
    glEnd();
  }
}

void setColor(int i) {
  switch (color_option) {

    case 0: {
      glColor4f(0.8f, 0.8f, 0.9f, 1.0f);
    } break;

    case 1: {
      colorRGBA col;
      colorTable.getRGB(Conf.grains[i].zncc1, &col);
      glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
    } break;

    case 2: {
      colorRGBA col;
      colorTable.getRGB(Conf.grains[i].zncc2, &col);
      glColor4f(col.r / 255.0, col.g / 255.0, col.b / 255.0, 1.0f);
    } break;

    default: {
      glColor4f(0.8f, 0.8f, 0.9f, 1.0f);
    } break;
  }
}

void drawParticles() {
  glLineWidth(1.0f);

  for (size_t i = 0; i < Conf.grains.size(); ++i) {
    double xc = Conf.grains[i].x;
    double yc = Conf.grains[i].y;
    double R = Conf.grains[i].R;

    setColor(i);
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
      double rot = Conf.grains[i].rot;
      glBegin(GL_LINES);
      glVertex2f(xc, yc);
      glVertex2f(xc + R * cos(rot), yc + R * sin(rot));
      glEnd();
    }
  }
}

void drawContacts() {}

void drawForces() {}

void drawDisplacements() {
  Conf.computeDisplacements(RefConf);
  
  glLineWidth(1.5f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  glBegin(GL_LINES);
  for (size_t i = 0; i < Conf.grains.size(); ++i) {
    double xc = Conf.grains[i].x;
    double yc = Conf.grains[i].y;
    double dx = Conf.grains[i].dx; // * vScale
    double dy = Conf.grains[i].dy; // * vScale
    
    glVertex2f(xc - dx, yc - dy);
    glVertex2f(xc, yc);
  }
  glEnd();
}

void drawFluctuations() {
  Conf.computeFluctuations(RefConf);
  
  glLineWidth(1.5f);
  glColor4f(1.0f, 0.0f, 0.0f, 1.0f);

  glBegin(GL_LINES);
  for (size_t i = 0; i < Conf.grains.size(); ++i) {
    double xc = Conf.grains[i].x;
    double yc = Conf.grains[i].y;
    double ddx = Conf.grains[i].ddx * vScale;
    double ddy = Conf.grains[i].ddy * vScale;
    
    glVertex2f(xc - ddx, yc - ddy);
    glVertex2f(xc, yc);
  }
  glEnd();
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

void try_to_readConf(int num, Configuration& CF, int& OKNum) {
  char file_name[256];
  int format = format_x_y_R_rot_NCC_NCC;  // format_x_y_R_rot
  sprintf(file_name, "CONF%04d", num);
  if (fileExists(file_name)) {
    std::cout << "Read " << file_name << std::endl;
    OKNum = num;
    CF.read(".", num, format);
  } else
    std::cout << file_name << " does not exist" << std::endl;
}

void menu(int num) {
  switch (num) {

    case 0:
      exit(0);
      break;
  };

  glutPostRedisplay();
}

void buildMenu() {
  int submenu1 = glutCreateMenu(menu);  // Particle Colors
  glutAddMenuEntry("None", 100);
  glutAddMenuEntry("Velocity Magnitude", 101);
  glutAddMenuEntry("Sum of Normal Contact Forces", 102);

  int submenu2 = glutCreateMenu(menu);  // Force Colors
  glutAddMenuEntry("None", 200);
  glutAddMenuEntry("Magnitude", 201);

  glutCreateMenu(menu);  // Main menu
  glutAddSubMenu("Particle Colors", submenu1);
  glutAddSubMenu("Force Colors", submenu2);
  glutAddSubMenu("Velocity Colors", submenu2);
  glutAddMenuEntry("Quit", 0);
}

// =====================================================================
// Main function
// =====================================================================

int main(int argc, char* argv[]) {

  if (argc == 1) {
    refConfNum = confNum = 1;
    ref_fixed = 1;
  } else if (argc == 2) {
    refConfNum = confNum = atoi(argv[1]);
    ref_fixed = 1;
  } else if (argc == 3) {
    refConfNum = atoi(argv[1]);
    confNum = atoi(argv[2]);
    ref_fixed = 0;
  }

  std::cout << "Reference Configuration: ";
  try_to_readConf(refConfNum, RefConf, refConfNum);
  std::cout << "Current Configuration: ";
  try_to_readConf(confNum, Conf, confNum);

  // ==== Init GLUT and create window
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA);
  glutInitWindowPosition(50, 50);
  glutInitWindowSize(width, height);
  main_window = glutCreateWindow("CONF VISUALIZER");

  // ==== Register callbacks
  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);

  // ==== Menu
  buildMenu();
  glutAttachMenu(GLUT_RIGHT_BUTTON);

  mouse_mode = NOTHING;

  glEnable(GL_BLEND);
  glBlendEquation(GL_FUNC_ADD);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // ==== Enter GLUT event processing cycle
  fit_view();
  glutMainLoop();
  return 0;
}
