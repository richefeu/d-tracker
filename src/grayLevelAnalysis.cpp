/************************************************************************************************/
/*                                     GRAY LEVEL ANALYSIS                                      */
/************************************************************************************************/

#include "moment.hpp"

void gray_level_analysis() {
  cout << "GRAY LEVEL ANALYSIS" << endl;
  read_image(im_index_ref, iref, true); // Read reference image

  std::vector<std::ofstream> file(num_grains);
  for (int igrain = 0; igrain < num_grains; igrain++) {
    char fname[256];
    sprintf(fname, "gray_level%d.txt", igrain);
    file[igrain].open(fname);
  }

  for (int num_image = ibeg; num_image <= iend; num_image += iinc) {
    read_image(im_index_current, num_image);

    for (int igrain = 0; igrain < num_grains; igrain++) {

      std::vector<double> val;
      double min = 65535.0 + 1.0;
      double max = 0.0;
      int xpixel, ypixel;
      for (size_t i = 0; i < grain[igrain].pattern.size(); i++) {
        xpixel = grain[igrain].pattern[i].dx + grain[igrain].refcoord_xpix;
        ypixel = grain[igrain].pattern[i].dy + grain[igrain].refcoord_ypix;
        double gray = (double)image[im_index_current][xpixel][ypixel];
        val.push_back(gray);
        if (min > gray)
          min = gray;
        if (max < gray)
          max = gray;
      }

      data_stat stat;
      moment(val, stat);

      file[igrain] << num_image << " " << stat.ave << " " << stat.adev << " " << stat.sdev << " " << stat.var << " "
                   << stat.skew << " " << stat.curt << " " << min << " " << max << '\n';
      // 1. image number
      // 2. average 
      // 3. absolut deviation
      // 4. standard deviation
      // 5. variance
      // 6. skewness
      // 7. curtosis
      // 8. min
      // 9. max
    }
  }
}