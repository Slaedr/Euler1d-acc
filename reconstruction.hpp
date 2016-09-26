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
	const std::vector<double>& x;						///< Cell centres
	const std::vector<std::vector<double>>& u;			///< Cell-centred variables
	std::vector<std::vector<double>>& dudx;				///< Cell-centred slopes
public:
	SlopeReconstruction(const int _N, const std::vector<double>& x, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx);
	virtual ~SlopeReconstruction();
	virtual void compute_slopes() = 0;
};

/// Just sets zero derivatives, used for 1st order scheme
class TrivialSlopeReconstruction : public SlopeReconstruction
{
public:
	TrivialSlopeReconstruction(const int _N, const std::vector<double>& x, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx);
	void compute_slopes();
};

/// Computes cell-centred derivatives using least-squares reconstruction
class LeastSquaresReconstruction : public SlopeReconstruction
{
public:
	LeastSquaresReconstruction(const int _N, const std::vector<double>& x, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx);
	void compute_slopes();
};


/// Interface for computing face values from cell-centred values
class FaceReconstruction
{
protected:
	const int N;										///< Number of cells
	const std::vector<double>& x;						///< Positions of cell centres
	const std::vector<std::vector<double>>& u;			///< Cell-centred variables
	const std::vector<std::vector<double>>& dudx;		///< Cell-centred slopes
	std::vector<std::vector<double>>& uleft;			///< Left value at each face
	std::vector<std::vector<double>>& uright;			///< Right value at each face
	std::string limiter;								///< String describing the limiter to use
public:
	FaceReconstruction(const int _N, const std::vector<double>& x, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, std::vector<std::vector<double>>& uleft,
			std::vector<std::vector<double>>& uright, std::string _limiter);
	virtual ~FaceReconstruction();
	virtual void compute_face_values() = 0;
};

class MUSCLReconstruction : public FaceReconstruction
{
	double k;											///< Controls order of reconstruction; people generally use 1/3
	const SlopeLimiter* lim;							///< Slope limiter to use
public:
	MUSCLReconstruction(const int _N, const std::vector<double>& x, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, std::vector<std::vector<double>>& uleft,
			std::vector<std::vector<double>>& uright, std::string _limiter, double _k);
	~MUSCLReconstruction();
	void compute_face_values();
};

#endif
