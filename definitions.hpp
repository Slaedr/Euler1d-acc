#ifndef __DEFITIONS_H
#define __DEFINITIONS_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#define ZERO_TOL 2.2e-16
#define NVARS 3

/// adiabatic index
const double g = 1.4;

/// Specific gas constant of air
//const double R = 287.1;
const double R = 1716;

const double Cv = R/(g-1.0);

void matprint(const std::vector<std::vector<double>>& mat);

#endif
