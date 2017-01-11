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
	#pragma acc routine seq
	virtual void compute_flux(double const *const uleft, double const *const uright, double *const flux) = 0;
	/// Compute flux from primitive variables at left and right of faces
	#pragma acc routine seq
	virtual void compute_flux_prim(double const *const uleft, double const *const uright, double *const flux) = 0;

	~InviscidFlux();
};

class LocalLaxFriedrichsFlux : public InviscidFlux
{
public:
	LocalLaxFriedrichsFlux();

	#pragma acc routine seq
	void compute_flux(double const *const uleft, double const *const uright, double *const flux);
	#pragma acc routine seq
	void compute_flux_prim(double const *const uleft, double const *const uright, double *const flux);
};

/// Van-Leer flux
class VanLeerFlux : public InviscidFlux
{
public:
	VanLeerFlux();

	#pragma acc routine seq
	void compute_flux(double const *const uleft, double const *const uright, double *const flux);
	#pragma acc routine seq
	void compute_flux_prim(double const *const uleft, double const *const uright, double *const flux);
};

#pragma acc routine seq
void compute_vanleerflux_prim(double const *const uleft, double const *const uright, double *const flux, const double g);

#endif

