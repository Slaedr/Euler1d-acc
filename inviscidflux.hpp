/** \file inviscidflux.hpp
 * \brief Approximate Riemann solvers for Euler equations
 * \author Aditya Kashi
 * \date 29 August 2016
 */

#ifndef __INVISCIDFLUX_H

#define __INVISCIDFLUX_H

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

/// All inviscid flux classes should inherit from this class
class InviscidFlux
{
public:
	virtual void compute_flux(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux) = 0;
};

class LocalLaxFriedrichsFlux : public InviscidFlux
{
public:
	void compute_flux(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux);
};

/// Van-Leer flux
class VanLeerFlux : public InviscidFlux
{
	std::vector<double> fluxL;
	std::vector<double> fluxR;
public:
	VanLeerFlux();
	void compute_flux(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux);
};

#endif

