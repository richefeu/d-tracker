#include "tracker.cpp"

int main(int argc, char **argv) {
  if (init(argc, argv) == 0)
    return 0;

  if (procedure == "correction_distortion")
    correction_distortion();
  else if (procedure == "correction_distortion_grid")
    correction_distortion_grid();
  else if (procedure == "find_subpixel_centers")
    find_subpixel_centers();
  else if (procedure == "pattern_quality")
    pattern_quality();
  else if (procedure == "pattern_fft")
    pattern_fft();
  else if (procedure == "post_process")
    post_process();
  else if (procedure == "visu_process")
    visu_process();
  else if (procedure == "gray_level_analysis")
    gray_level_analysis();
  else
    particle_tracking(); // The default procedure
  return 0;
}