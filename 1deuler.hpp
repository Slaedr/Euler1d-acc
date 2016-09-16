/** \file 1deuler.hpp
 * \brief Code to control the solution process for time-accurate 1D Euler equations.
 * \author Aditya Kashi
 * \date 29 August 2016
 */

#ifndef __1DEULER_H
#define __1DEULER_H

#ifndef __INVISCIDFLUX_H
#include "inviscidflux.hpp"
#endif

/// Base class for Euler solution processes
class Euler1d
{
protected:
	int N;										///< Number of cells
	std::vector<double> x;						///< Cell centers
	std::vector<double> dx;						///< (1D) Size of each cell
	std::vector<double> vol;					///< (3D) Volume of each cell
	std::vector<double> nodes;					///< Mesh nodes
	double domlen;								///< Physical length of the domain
	std::vector<double> A;						///< Cross-sectional areas at cell centers
	std::vector<std::vector<double>> u;			///< Unknowns - u[i][0] is density of cell i and so on
	std::vector<std::vector<double>> uleft;		///< Left state of each face
	std::vector<std::vector<double>> uright;	///< Right state at each face
	std::vector<std::vector<double>> dudx;		///< Slope of conserved variables in each cell
	std::vector<std::vector<double>> res;		///< residual
	int bcL;									///< left BC type
	int bcR;									///< right BC type
	std::vector<double> bcvalL;					///< left boundary value
	std::vector<double> bcvalR;					///< right boundary value
	InviscidFlux* flux;							///< Inviscid flux computation context
	SlopeReconstruction* cslope;				///< Slope reconstruction context
	FaceReconstruction* rec;					///< Context responsible for computation of face values of flow variables from their cell-centred values
	double cfl;									///< CFL number			

public:
	Euler1d(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs, std::string inviscid_flux, double CFL);

	virtual ~Euler1d();

	/// Generates a grid depending on 
	/** \param type If type == 0, a uniform grid is generated and used. If type == 1, then grid points need to passed to the function in
	 * \param pointlist an array of positions of mesh points.
	 *
	 * Note that the domain is assumed to start at x=0.
	 */
	void generate_mesh(int type, const std::vector<double>& pointlist);

	/// Set cross-sectional areas
	void set_area(int type, std::vector<double>& cellCenteredAreas);

	void compute_slopes();

	void compute_face_values();

	void compute_inviscid_fluxes();

	void compute_source_term();

	void apply_boundary_conditions();

	void postprocess(std::string outfilename);
};

/// Explicit RK solver for time-dependent 1D Euler equations
class Euler1dExplicit : public Euler1d
{
	double ftime;								///< Physical time for which to simulate
	std::vector<double> maxWaveSpeed;			///< for computing time steps
	int temporalOrder;							///< desired temporal order of accuracy

public:
	Euler1dExplicit(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs, std::string inviscidFlux, double cfl, 
			double fTime, int temporal_order);

	void run();
};

/// Explicit RK solver for steady-state 1D Euler
class Euler1dSteadyExplicit : public Euler1d
{
	double tol;
	int maxiter;
	std::vector<double> maxWaveSpeed;			///< for computing time steps

public:
	Euler1dSteadyExplicit(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs, std::string inviscidFlux, double cfl, 
			double toler, int max_iter);

	void run();
};

#endif
