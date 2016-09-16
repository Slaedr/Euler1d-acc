/** \file reconstruction.hpp
 * \brief Schemes for computing the derivatives by reconstructio, and for computing face values from cell-centred values.
 * \author Aditya Kashi
 * \date 16 September 2016
 */

#ifndef __RECONSTRUCTION_H

#ifndef __DEFINITIONS_H
#include "definitions.hpp"
#endif

/// Interface for slope reconstruction classes
class SlopeReconstruction
{
protected:
	const int N;										///< Number of cells
	const std::vector<std::vector<double>>& u;			///< Cell-centred variables
	std::vector<std::vector<double>>& dudx;				///< Cell-centred slopes
public:
	SlopeReconstruction(const int _N, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx);
	virtual ~SlopeReconstruction();
	virtual void compute_slopes() = 0;
};

/// Just sets zero derivatives, used for 1st order scheme
class TrivialSlopeReconstruction : public SlopeReconstruction
{
public:
	void compute_slopes();
};

/// Interface for computing face values from cell-centred values
class FaceReconstruction
{
protected:
	const int N;										///< Number of cells
	const std::vector<std::vector<double>>& u;			///< Cell-centred variables
	const std::vector<std::vector<double>>& dudx;		///< Cell-centred slopes
	std::vector<std::vector<double>>& uleft;			///< Left value at each face
	std::vector<std::vector<double>>& uright;			///< Right value at each face
public:
	FaceReconstruction(const int _N, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, std::vector<std::vector<double>>& uleft,
			std::vector<std::vector<double>>& uright);
	virtual ~FaceReconstruction();
	virtual void compute_face_values() = 0;
};

#endif
