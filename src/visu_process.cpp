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
  cout << "VISU PROCESS" << endl;

  int num_image;
  char fname[256];
  vector<grain_type_2D> grain_ref;

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

    thumbnail thumb(image, 0, image_div);
    ColorTable ct;
    ct.setMinMax(colorMin, colorMax);
    if (num_image == ibeg) {
      ct.SavePpm("ColorMap.ppm");
    }
    colorRGBA cRGBA;

    if (visu_mode == "continuum") { /// ============================================================================
      cout << "CONTINUUM" << endl;
      vector<Point2D> points;
      for (int i = 0; i < num_grains; ++i) {
        double x_centre = grain[i].refcoord_xpix + grain[i].dx;
        double y_centre = grain[i].refcoord_ypix + grain[i].dy;

        points.push_back(Point2D(nearest(x_centre), nearest(y_centre)));
      }

      // MESH
      vector<triangle_def> triangles;
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
      vector<double> data;

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

        cout << "exx_mean = " << exx_mean << endl;
        cout << "eyy_mean = " << eyy_mean << endl;
        cout << "exy_mean = " << exy_mean << endl;

        vector<point_def> point_defs;
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
        cerr << "nx_grid_visu not compatible with num_grains" << endl;
        return;
      }
      if (nx_grid_visu < npts) {
        cerr << "nx_grid_visu too small" << endl;
        return;
      }
      if (ny_grid_visu < npts) {
        cerr << "ny_grid_visu too small" << endl;
        return;
      }

      vector<double> Fxx(num_grains), Fxy(num_grains), Fyx(num_grains), Fyy(num_grains);
      vector<double> Exx(num_grains), Exy(num_grains), Eyx(num_grains), Eyy(num_grains);
      vector<double> XData(npts), YData(npts), UxData(npts), UyData(npts);
      vector<double> data;

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
      }      // End strain computation
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
        LOG("ColorScale / " << visu_output << " / im = " << num_image << " / min = " << colorMin
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
      cout << "COLOR_DISKS" << endl;
      RGB col;
      vector<double> data;
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
        SHOW(colorMin);
        SHOW(colorMax);
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
    //thumb.writePpm(fname, visu_alpha);
    thumb.writeTiff(fname, visu_alpha);
  } // End loop over files (image and dic_out_x.txt)
}
