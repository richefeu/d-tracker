#ifndef SHELLFILES_HPP
#define SHELLFILES_HPP

#include <iostream>
#include <fstream>
#include <vector>

#include "Data.hpp"

void saveShells(const char *name, std::vector<Shell> & shells) {
  std::ofstream outFile(name);
  outFile << shells.size() << std::endl;
  for (size_t i = 0; i < shells.size() ; i++) {
    outFile << shells[i].x << " " << shells[i].y << " " << shells[i].z << " "
      << shells[i].roty << " " << shells[i].rotz << " " 
      << shells[i].r_in << " " << shells[i].r_out << " " << shells[i].h << " " 
      << shells[i].best_match << std::endl;
  }
}

void readShells(const char *name, std::vector<Shell> & shells) {
  std::ifstream inFile(name);
  size_t nb = 0;
  inFile >> nb;
  shells.resize(nb);
  for (size_t i = 0; i < shells.size() ; i++) {
    inFile >> shells[i].x >> shells[i].y >> shells[i].z 
      >> shells[i].roty >> shells[i].rotz  
      >> shells[i].r_in >> shells[i].r_out >> shells[i].h >> shells[i].best_match;
  }
}

#endif /* end of include guard: SHELLFILES_HPP */
