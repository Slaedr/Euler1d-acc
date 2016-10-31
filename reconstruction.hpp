/** \file reconstruction.hpp
 * \brief Schemes for computing the derivatives by reconstructio, and for computing face values from cell-centred values.
 * \author Aditya Kashi
 * \date 16 September 2016
 */

#ifndef __RECONSTRUCTION_H

#ifndef __DEFINITIONS_H
#include "definitions.hpp"
#endif

#ifndef __LIMITERS_H
#include "limiters.hpp"
#endif

/// Interface for slope reconstruction classes
class SlopeReconstruction
{
protected:
	const int N;										///< Number of cells
	double const *const x;						///< Cell centres
	double const *const dx;						///< Cell widths
	double const *const *const u;			///< Cell-centred variables
	double * const * const dudx;				///< Cell-centred slopes
public:
	SlopeReconstruction(const int _N, double const *const x, double const *const dx, double const *const *const _u, double * const * const _dudx);
	virtual ~SlopeReconstruction();
	virtual void compute_slopes() = 0;
};

/// Just sets zero derivatives, used for 1st order scheme
class TrivialSlopeReconstruction : public SlopeReconstruction
{
public:
	TrivialSlopeReconstruction(const int _N, double const *const x, double const *const dx, double const *const *const _u, double * const * const _dudx);
	void compute_slopes();
};

/// Computes cell-centred derivatives using least-squares reconstruction
class LeastSquaresReconstruction : public SlopeReconstruction
{
public:
	LeastSquaresReconstruction(const int _N, double const *const x, double const *const dx, double const *const *const _u, double * const * const _dudx);
	void compute_slopes();
};

/// Computes cell-centred derivatives using central difference
class CentralDifferenceReconstruction : public SlopeReconstruction
{
public:
	CentralDifferenceReconstruction(const int _N, double const *const x, double const *const dx, double const *const *const _u, double * const * const _dudx);
	void compute_slopes();
};

/// Computes TVD limited slopes
class TVDSlopeReconstruction : public SlopeReconstruction
{
	const SlopeLimiter1* lim;				///< Slope limiter to use
public:
	TVDSlopeReconstruction(const int _N, double const *const x, double const *const dx, double const *const *const _u, 
			double * const * const _dudx, std::string limiter);
	~TVDSlopeReconstruction();
	void compute_slopes();
};


/// Interface for computing face values from cell-centred values
class FaceReconstruction
{
protected:
	const int N;										///< Number of cells
	double const *const x;						///< Positions of cell centres
	double const *const *const u;			///< Cell-centred variables
	double const *const *const dudx;		///< Cell-centred slopes
	double * const * const uleft;			///< Left value at each face
	double * const * const uright;			///< Right value at each face
public:
	FaceReconstruction(const int _N, double const *const x, double const *const *const _u, double const *const *const _dudx, double * const * const uleft,
			double * const * const uright);
	virtual ~FaceReconstruction();
	virtual void compute_face_values() = 0;
};

class MUSCLReconstruction : public FaceReconstruction
{
	double k;											///< Controls order of reconstruction; people generally use 1/3
	std::string limiter;								///< String describing the limiter to use
	const Limiter* lim;									///< Slope limiter to use
public:
	MUSCLReconstruction(const int _N, double const *const x, double const *const *const _u, double const *const *const _dudx, double * const * const uleft,
			double * const * const uright, std::string _limiter, double _k);
	~MUSCLReconstruction();
	void compute_face_values();
};

class MUSCLReconstructionG : public FaceReconstruction
{
	double k;											///< Controls order of reconstruction; people generally use 1/3
	std::string limiter;								///< String describing the limiter to use
	const SlopeLimiter1* lim;							///< Slope limiter to use
public:
	MUSCLReconstructionG(const int _N, double const *const x, double const *const *const _u, double const *const *const _dudx, double * const * const uleft,
			double * const * const uright, std::string _limiter, double _k);
	~MUSCLReconstructionG();
	void compute_face_values();
};

class LinearReconstruction : public FaceReconstruction
{
public:
	LinearReconstruction(const int _N, double const *const x, double const *const *const _u, double const *const *const _dudx, double * const * const uleft,
			double * const * const uright);
	void compute_face_values();
};

/*class PureMUSCLReconstruction : public FaceReconstruction
{
	double k;											///< Controls order of reconstruction; people generally use 1/3
	const SlopeLimiter* lim;							///< Slope limiter to use
public:
	MUSCLReconstruction(const int _N, double const *const x, double const *const *const _u, double const *const *const _dudx, double * const * const uleft,
			double * const * const uright, std::string _limiter, double _k);
	~MUSCLReconstruction();
	void compute_face_values();
};*/

#endif
