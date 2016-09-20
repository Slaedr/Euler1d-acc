/** \file reconstruction.cpp
 * \brief Implementation of reconstruction schemes.
 * \author Aditya Kashi
 * \date
 */

#include "reconstruction.hpp"

SlopeReconstruction::SlopeReconstruction(const int _N, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: N(_N), u(_u), dudx(_dudx)
{ }
	
SlopeReconstruction::~SlopeReconstruction()
{ }

TrivialSlopeReconstruction::TrivialSlopeReconstruction(const int _N, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: SlopeReconstruction(_N, _u, _dudx)
{ }

void TrivialSlopeReconstruction::compute_slopes()
{
	for(int i = 1; i <= N; i++)
		for(int j = 0; j < NVARS; j++)
			dudx[i][j] = 0;
}
	
FaceReconstruction(const int _N, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, std::vector<std::vector<double>>& _uleft,
			std::vector<std::vector<double>>& _uright, std::string _limiter)
	: N(_N), u(_u), dudx(_dudx), uleft(_uleft), uright(_uright), limiter(_limiter)
{ }
	
FaceReconstruction::~FaceReconstruction()
{ }

MUSCLReconstruction(const int _N, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, std::vector<std::vector<double>>& uleft, 
		std::vector<std::vector<double>>& uright, std::string _limiter, double _k)
	: FaceReconstruction(_N, _u, _dudx, uleft, uright, _limiter)
{
	if(limiter == "vanalbada")
		lim = new VanAlbadaLimiter(SMALL_NUMBER);
	else
		lim = new VanAlbadaLimiter(SMALL_NUMBER);
}

MUSCLReconstruction::~MUSCLReconstruction()
{
	delete lim;
}
