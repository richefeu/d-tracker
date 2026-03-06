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
  
  cout << "Number of points for equiprojectivity computations: " << List_Grains_i.size() << '\n';
}

// La fonction a minimiser si on veut utiliser le critère d'équiprojectivité.
// Pour il faut que les déplacements entre les image correspondent réellement à un mouvement
// de solide rigide (déplacement d'une grande plaque).
//
// Remarque suite aux tests de Gael :
// La minimisation de sum(|err|) se comporte mal! (solution très sensible aux perturbations de départ)
// Il vaut mieux minimiser sum(err). Mais on préfère finalement minimiser max(|err|)
double disto_to_minimize_Equiproj(vector<double> &X) {
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

      //undistor(&X[0], (double)(grain[igrain].refcoord_xpix), (double)(grain[igrain].refcoord_ypix), A1x, A1y);
      A1x = (double)(grain[igrain].refcoord_xpix);
      A1y = (double)(grain[igrain].refcoord_ypix);
      undistor(&X[0], (double)(grain[igrain].refcoord_xpix + dx_corrDisto[i_image][igrain]),
               (double)(grain[igrain].refcoord_ypix + dy_corrDisto[i_image][igrain]), A2x, A2y);
      //undistor(&X[0], (double)(grain[jgrain].refcoord_xpix), (double)(grain[jgrain].refcoord_ypix), B1x, B1y);
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
      val = (ABx * (BBx - AAx)) + (ABy * (BBy - AAy)) ; 
      func += val * val;
      //if (val > func) func = val;
    }
  }
   // return fabs(func);
   return (func);
}

double Equiproj_Error_Computation(vector<double> &X) {
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
  cout << "CORRECTION OF DISTORTION" << endl;

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

  ofstream logfile("corrDisto.log");

  read_image(0, image_numbers_corrDisto[0], true); // Read reference image
  iref = image_numbers_corrDisto[0];
  do_precomputations();
  precompute_paires();

  if (make_images)
    create_image(0);

  cout << "Tracking displacements... " << endl << flush;
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
      { loadbar(++progress, grain.size()); }
    }
    cerr << endl;
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
      { loadbar(++progress, grain.size()); }
    }
    cerr << endl;
    fprintf(stdout, "[DONE in %f seconds]\n", get_time() - tbeg);

    int nbFailed = 0;
    for (igrain = 0; igrain < num_grains; igrain++) {
      dx_corrDisto[i_image][igrain] = grain[igrain].dx;
      dy_corrDisto[i_image][igrain] = grain[igrain].dy;
      NCC_subpix_corrDisto[i_image][igrain] = grain[igrain].NCC_subpix;
      if (grain[igrain].NCC_subpix < NCC_min)
        nbFailed++;
    }
    cout << "Failure rate = " << 100.0 * (double)nbFailed / (double)num_grains << "%" << endl;

    // save_grains(i_image);
    save_grains(image_numbers_corrDisto[i_image]);

  } // for loop i_image
  cout << "done." << endl;

  // The parameters (initial guess) and perturbations are set by the user
  if (disto_parameters[0] < 0.0)
    disto_parameters[0] = 0.5 * (double)dimx;
  if (disto_parameters[1] < 0.0)
    disto_parameters[1] = 0.5 * (double)dimy;

  index_ini = Equiproj_Error_Computation(disto_parameters);

  // Minimization
  int npairs = List_Grains_i.size();
  cout << "Numbers of pairs used to undistor : " << npairs << endl;

  cout << "Minimizing equiprojectivity criterion... " << flush;
  Powell<double(vector<double> &)> powell(disto_to_minimize_Equiproj, 1e-12);
  disto_parameters = powell.minimize(disto_parameters, disto_parameters_perturb);
  cout << "done." << endl << endl;

  index_end = Equiproj_Error_Computation(disto_parameters);
  index_gain = index_ini / index_end - 1.0;
  logfile << scientific << setprecision(10);
  logfile << "fake_undistor            " << fake_undistor << endl;
  logfile << "Initial distorsion index " << index_ini << endl;
  logfile << "Final distorsion index   " << index_end << endl;
  logfile << "Index gain (>0 is good)  " << index_gain << endl;
  logfile << "xc_corrDistor            " << disto_parameters[0] << endl;
  logfile << "yc_corrDistor            " << disto_parameters[1] << endl;
  logfile << "K1_corrDistor            " << disto_parameters[2] << endl;
  logfile << "K2_corrDistor            " << disto_parameters[3] << endl;
  logfile << "K3_corrDistor            " << disto_parameters[4] << endl;
  logfile << "P1_corrDistor            " << disto_parameters[5] << endl;
  logfile << "P2_corrDistor            " << disto_parameters[6] << endl;
  logfile << "P3_corrDistor            " << disto_parameters[7] << endl;

  cout << scientific << setprecision(10);
  cout << "fake_undistor            " << fake_undistor << endl;
  cout << "Initial distorsion index " << index_ini << endl;
  cout << "Final distorsion index   " << index_end << endl;
  cout << "Index gain (>0 is good)  " << index_gain << endl;
  cout << "xc_corrDistor            " << disto_parameters[0] << endl;
  cout << "yc_corrDistor            " << disto_parameters[1] << endl;
  cout << "K1_corrDistor            " << disto_parameters[2] << endl;
  cout << "K2_corrDistor            " << disto_parameters[3] << endl;
  cout << "K3_corrDistor            " << disto_parameters[4] << endl;
  cout << "P1_corrDistor            " << disto_parameters[5] << endl;
  cout << "P2_corrDistor            " << disto_parameters[6] << endl;
  cout << "P3_corrDistor            " << disto_parameters[7] << endl;

  ofstream errorFile("errors.dat");
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
      if (fabs(xerr) > xerrormax) xerrormax = fabs(xerr);
      if (fabs(yerr) > yerrormax) yerrormax = fabs(yerr);
      errorFile << xd << " " << yd << " " << xerr << " " << yerr << endl;
    }
    errorFile << endl;
  }
  logfile << "xerrormax                " << xerrormax << endl;
  logfile << "yerrormax                " << yerrormax << endl;
  cout << "xerrormax                " << xerrormax << endl;
  cout << "yerrormax                " << yerrormax << endl;
  
  ofstream errorBoxFile("errorbox.dat");
  yd = 0;
  for (xd = 0; xd < dimx; xd += dxd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
  }
  xd = dimx - 1;
  for (yd = 0; yd < dimy; yd += dyd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
  }
  yd = dimy - 1;
  for (xd = dimx - 1; xd >= 0; xd -= dxd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
  }
  xd = 0;
  for (yd = dimy - 1; yd >= 0; yd -= dyd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
  }

  xd = yd = 0;
  undistor(&disto_parameters[0], xd, yd, xu, yu);
  errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
}

void correction_distortion_grid() {
  cout << "CORRECTION OF DISTORSION -- GRID METHOD" << endl;
  read_image(im_index_ref, grid_image_name.c_str(), true); // Read grid image (to get dimx and dimy)
  ofstream logfile("corrDisto.log");
  double index = disto_to_minimize_grid(disto_parameters);
  logfile << "Initial distorsion index " << index << endl;
  cout << "Initial distorsion index " << index << endl;

  // The parameters (initial guess) and perturbations are set by the user
  if (disto_parameters[0] < 0.0)
    disto_parameters[0] = 0.5 * (double)dimx;
  if (disto_parameters[1] < 0.0)
    disto_parameters[1] = 0.5 * (double)dimy;

  // Minimization
  Powell<double(vector<double> &)> powell(disto_to_minimize_grid, 1e-8);
  disto_parameters = powell.minimize(disto_parameters, disto_parameters_perturb);

  index = disto_to_minimize_grid(disto_parameters);
  logfile << "Final distorsion index   " << index << endl;
  logfile << "xc_corrDistor            " << disto_parameters[0] << endl;
  logfile << "yc_corrDistor            " << disto_parameters[1] << endl;
  logfile << "K1_corrDistor            " << disto_parameters[2] << endl;
  logfile << "K2_corrDistor            " << disto_parameters[3] << endl;
  logfile << "K3_corrDistor            " << disto_parameters[4] << endl;
  logfile << "P1_corrDistor            " << disto_parameters[5] << endl;
  logfile << "P2_corrDistor            " << disto_parameters[6] << endl;
  logfile << "P3_corrDistor            " << disto_parameters[7] << endl;

  cout << "Final distorsion index   " << index << endl;
  cout << "xc_corrDistor            " << disto_parameters[0] << endl;
  cout << "yc_corrDistor            " << disto_parameters[1] << endl;
  cout << "K1_corrDistor            " << disto_parameters[2] << endl;
  cout << "K2_corrDistor            " << disto_parameters[3] << endl;
  cout << "K3_corrDistor            " << disto_parameters[4] << endl;
  cout << "P1_corrDistor            " << disto_parameters[5] << endl;
  cout << "P2_corrDistor            " << disto_parameters[6] << endl;
  cout << "P3_corrDistor            " << disto_parameters[7] << endl;

  ofstream gridFile("grid.dat");
  double xd, yd, xu, yu;
  for (int iy = 0; iy < ny_grid_disto; iy++) {
    for (int ix = 0; ix < nx_grid_disto; ix++) {
      int igrain = iy * nx_grid_disto + ix;
      xd = grain[igrain].refcoord_xpix;
      yd = grain[igrain].refcoord_ypix;
      undistor(&disto_parameters[0], xd, yd, xu, yu);
      gridFile << xd << " " << yd << " " << xu << " " << yu << endl;
    }
    gridFile << endl;
  }
  gridFile << endl;
  for (int ix = 0; ix < nx_grid_disto; ix++) {
    for (int iy = 0; iy < ny_grid_disto; iy++) {
      int igrain = iy * nx_grid_disto + ix;
      xd = grain[igrain].refcoord_xpix;
      yd = grain[igrain].refcoord_ypix;
      undistor(&disto_parameters[0], xd, yd, xu, yu);
      gridFile << xd << " " << yd << " " << xu << " " << yu << endl;
    }
    gridFile << endl;
  }

  ofstream errorFile("errors.dat");
  double dxd = (double)dimx / 100.0;
  double dyd = (double)dimy / 100.0;
  if (dxd < 1.0)
    dxd = 1.0;
  if (dyd < 1.0)
    dyd = 1.0;
  for (xd = 0; xd < dimx; xd += dxd) {
    for (yd = 0; yd < dimy; yd += dyd) {
      undistor(&disto_parameters[0], xd, yd, xu, yu);
      errorFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
    }
    errorFile << endl;
  }

  ofstream errorBoxFile("errorbox.dat");
  yd = 0;
  for (xd = 0; xd < dimx; xd += dxd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
  }
  xd = dimx - 1;
  for (yd = 0; yd < dimy; yd += dyd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
  }
  yd = dimy - 1;
  for (xd = dimx - 1; xd >= 0; xd -= dxd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
  }
  xd = 0;
  for (yd = dimy - 1; yd >= 0; yd -= dyd) {
    undistor(&disto_parameters[0], xd, yd, xu, yu);
    errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;
  }

  xd = yd = 0;
  undistor(&disto_parameters[0], xd, yd, xu, yu);
  errorBoxFile << xd << " " << yd << " " << xu - xd << " " << yu - yd << endl;

  // show_distor();
  undistor_image(grid_image_name.c_str(), "undisto.png");
}

// To force horizontal lines to be horizontal, and vertical lines to be vertical
double disto_to_minimize_grid(vector<double> &X) {
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

double SolidMotionError_to_minimize(vector<double> &X) {
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
  vector<double> X(5);
  X[0] = angle;
  X[1] = xc;
  X[2] = yc;
  X[3] = xtrans;
  X[4] = ytrans;

  vector<double> dX(5);
  dX[0] = 0.0; // 1.0e-6;
  dX[1] = 0.0; // 1.0e-6;
  dX[2] = 0.0; // 1.0e-6;
  dX[3] = 1.0e-6;
  dX[4] = 1.0e-6;

  // Minimization
  Powell<double(vector<double> &)> powell(SolidMotionError_to_minimize, 1e-8);
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
