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

  vector<grain_type_2D> gr;
  vector<sensor_data> sensor;

  vector<processed_point> pt;      // The points corresponding to sample elements
  vector<processed_point> corners; // 4 corners of the loading parallepiped
  vector<processed_point> fix_pt0; // 2 first fixe points on the first dic_out file
  vector<processed_point> fix_pt;  // All fixe points

  process_data() : scaleFactor(1.0) {}
};

void read_sensor_data(const char *name, process_data &data) {
  data.sensor.clear();
  ifstream file(name);
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
  ifstream file(name);
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
  ofstream conf(name);

  conf << data.pt.size() << endl;
  conf << scientific << setprecision(5);
  for (size_t i = 0; i < data.pt.size(); i++) {
    conf << data.pt[i].x << ' ' << data.pt[i].y << ' ' << data.pt[i].rot << ' ' << data.pt[i].radius << endl;
  }

  conf << "#CORNERS  The four corners" << endl;
  conf << data.corners.size() << endl;
  for (size_t i = 0; i < data.corners.size(); i++) {
    conf << data.corners[i].x << ' ' << data.corners[i].y << ' ' << data.corners[i].rot << endl;
  }

  conf << "#FIX_PTS  The fix points" << endl;
  conf << data.fix_pt.size() << endl;
  for (size_t i = 0; i < data.fix_pt.size(); i++) {
    conf << data.fix_pt[i].x << ' ' << data.fix_pt[i].y << ' ' << data.fix_pt[i].rot << endl;
  }

  conf << "#SCALES_  The scale factor expressed in length per pixel (e.g. m/pix)" << endl;
  conf << data.scaleFactor << endl;

  conf << "#FILES__  The two files used to correlate" << endl;
  conf << data.imFile_FROM << ' ' << data.imFile_TO << endl;
}

void post_process() {
  cout << "POST PROCESS" << endl;
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
      cerr << "Cannot rescale, at least 2 fixed points need to be set" << endl;
      exit(0);
    }
    if (scaling_distance == 0.0) {
      cerr << "Cannot rescale, scaling_distance need to be set" << endl;
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
