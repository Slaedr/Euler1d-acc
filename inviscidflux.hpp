/** \file inviscidflux.hpp
 * \brief Approximate Riemann solvers for Euler equations
 * \author Aditya Kashi
 * \date 29 August 2016
 */

#ifndef __INVISCIDFLUX_H
#define __INVISCIDFLUX_H

#ifndef __DEFINITIONS_H
#include "definitions.hpp"
#endif

/// All inviscid flux classes should inherit from this class
class InviscidFlux
{
public:
	/// Compute flux from converved variables at left and right of faces
	virtual void compute_flux(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux) = 0;
	/// Compute flux from primitive variables at left and right of faces
	virtual void compute_flux_prim(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux) = 0;
};

class LocalLaxFriedrichsFlux : public InviscidFlux
{
public:
	void compute_flux(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux);
	void compute_flux_prim(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux);
};

/// Van-Leer flux
class VanLeerFlux : public InviscidFlux
{
	std::vector<double> fluxL;
	std::vector<double> fluxR;
public:
	VanLeerFlux();
	void compute_flux(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux);
	void compute_flux_prim(const std::vector<double>& uleft, const std::vector<double>& uright, std::vector<double>& flux);
};

#endif

