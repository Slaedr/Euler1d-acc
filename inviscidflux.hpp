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
	InviscidFlux();

	/// Compute flux from converved variables at left and right of faces
	virtual void compute_flux(double const *const uleft, double const *const uright, double const *const flux) = 0;
	/// Compute flux from primitive variables at left and right of faces
	virtual void compute_flux_prim(double const *const uleft, double const *const uright, double const *const flux) = 0;

	~InviscidFlux();
};

class LocalLaxFriedrichsFlux : public InviscidFlux
{
public:
	LocalLaxFriedrichsFlux();

	void compute_flux(double const *const uleft, double const *const uright, double const *const flux);
	void compute_flux_prim(double const *const uleft, double const *const uright, double const *const flux);
};

/// Van-Leer flux
class VanLeerFlux : public InviscidFlux
{
public:
	VanLeerFlux();

	void compute_flux(double const *const uleft, double const *const uright, double const *const flux);
	void compute_flux_prim(double const *const uleft, double const *const uright, double const *const flux);
};

#endif

