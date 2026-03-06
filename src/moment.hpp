#ifndef MOMENT_HPP_AAEE01E9
#define MOMENT_HPP_AAEE01E9

#include <cmath>
#include <vector>

struct data_stat {
  double ave;
  double adev;
  double sdev;
  double var;
  double skew;
  double curt;
};

void moment(std::vector<double> &data, data_stat &stat) {
  size_t n = data.size();

  if (n <= 1)
    std::cerr << "@moment, number of data too small (min = 2)\n";
  double s = 0.0;
  for (size_t j = 0; j < n; j++)
    s += data[j];
  stat.ave = s / (double)n;
  stat.adev = stat.var = stat.skew = stat.curt = 0.0;
  double ep = 0.0;
  double p = 0.0;
  for (size_t j = 0; j < n; j++) {
    stat.adev += fabs(s = data[j] - stat.ave);
    ep += s;
    stat.var += (p = s * s);
    stat.skew += (p *= s);
    stat.curt += (p *= s);
  }
  stat.adev /= (double)n;
  stat.var = (stat.var - ep * ep / (double)n) / (double)(n - 1);
  stat.sdev = sqrt(stat.var);
  if (stat.var != 0.0) {
    stat.skew /= (double)(n * stat.var * stat.sdev);
    stat.curt = stat.curt / (double)(n * stat.var * stat.var) - 3.0;
  } else
    std::cerr << "@moment, No skew/kurtosis when variance = 0 (in moment)\n";
}

#endif /* end of include guard: MOMENT_HPP_AAEE01E9 */
