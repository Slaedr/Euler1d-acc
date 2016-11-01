#ifndef __DEFITIONS_H
#define __DEFINITIONS_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

#define ZERO_TOL 2.2e-16
#define SMALL_NUMBER 1e-10

#define NVARS 3

const double PI = 3.14159265359;

/// adiabatic index
const double g = 1.4;

/// Specific gas constant of air
//const double R = 287.1;
const double R = 1716;

const double Cv = R/(g-1.0);

/// Constant for MUSCL reconstruction
const double muscl_k = 1.0/3.0;

void matprint(const std::vector<std::vector<double>>& mat);

// OpenACC constants
#define NVIDIA_VECTOR_LENGTH 32

#endif
