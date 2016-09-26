/** \file reconstruction.cpp
 * \brief Implementation of reconstruction schemes.
 * \author Aditya Kashi
 * \date
 */

#include "reconstruction.hpp"

SlopeReconstruction::SlopeReconstruction(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: N(_N), x(_x), u(_u), dudx(_dudx)
{ }
	
SlopeReconstruction::~SlopeReconstruction()
{ }

TrivialSlopeReconstruction::TrivialSlopeReconstruction(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: SlopeReconstruction(_N, _x, _u, _dudx)
{ }

void TrivialSlopeReconstruction::compute_slopes()
{
	for(int i = 1; i <= N; i++)
		for(int j = 0; j < NVARS; j++)
			dudx[i][j] = 0;
}

LeastSquaresReconstruction::LeastSquaresReconstruction(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: SlopeReconstruction(_N, _x, _u, _dudx)
{ }

void LeastSquaresReconstruction::compute_slopes()
{
	// implement least-squares reconstruction
	double num, denom;
	for(int i = 1; i <= N; i++)
		for(int j = 0; j < NVARS; j++)
		{
			num = (u[i-1][j]-u[i][j])*(x[i-1]-x[i]) + u[i+1][j]-u[i][j]*(x[i+1]-x[i]);
			denom = (x[i-1]-x[i])*(x[i-1]-x[i]) + (x[i+1]-x[i])*(x[i+1]-x[i]);
			dudx[i][j] = num/denom;
		}
}
	
FaceReconstruction::FaceReconstruction(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, 
		std::vector<std::vector<double>>& _uleft, std::vector<std::vector<double>>& _uright, std::string _limiter)
	: N(_N), x(_x), u(_u), dudx(_dudx), uleft(_uleft), uright(_uright), limiter(_limiter)
{ }
	
FaceReconstruction::~FaceReconstruction()
{ }

MUSCLReconstruction::MUSCLReconstruction(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, 
		std::vector<std::vector<double>>& uleft, std::vector<std::vector<double>>& uright, std::string _limiter, double _k)
	: FaceReconstruction(_N, _x, _u, _dudx, uleft, uright, _limiter)
{
	if(limiter == "vanalbada")
		lim = new VanAlbadaLimiter(SMALL_NUMBER);
	else
	{
		lim = new NoLimiter(SMALL_NUMBER);
		std::cout << "MUSCLReconstruction: Caution: not using any limiter.\n";
	}
}

MUSCLReconstruction::~MUSCLReconstruction()
{
	delete lim;
}

void MUSCLReconstruction::compute_face_values()
{
	// iterate over interior face
	int i, j, k;
	double sminus, splus, delminus, delplus;
	for(i = 1; i < N; i++)
	{
		// get deltas
		for(j = 0; j < NVARS; j++)
		{
			delminus = 2.0*dudx[i][j]*(x[i+1]-x[i]);
			delplus = 2.0*dudx[i+1][j]*(x[i+1]-x[i]);

			sminus = lim->limiter_function(delminus, u[i+1][j]-u[i][j]);
			splus = lim->limiter_function(delplus, u[i+1][j]-u[i][j]);

			uleft[i][j] = u[i][j] + sminus/4.0*( (1-k*sminus)*delminus + (1+k*sminus)*(u[i+1][j]-u[i][j]) );
			uright[i][j] = u[i+1][j] - splus/4.0*( (1-k*splus)*delplus + (1+k*splus)*(u[i+1][j]-u[i][j]) );

			if(fabs(dudx[i][j]) < 1e-15)
				uleft[i][j] = u[i][j];
			if(fabs(dudx[i+1][j]) < 1e-15)
				uright[i][j] = u[i+1][j];
		}
	}

	// get uright at 0 and uleft at N (boundary faces)
	for(j = 0; j < NVARS; j++)
	{
		delplus = 2*dudx[1][j]*(x[1]-x[0]);
		splus = lim->limiter_function(delplus, u[1][j]-u[0][j]);
		uright[0][j] = u[1][j] - splus/4.0*( (1-k*splus)*delplus + (1+k*splus)*(u[1][j]-u[0][j]) );

		if(fabs(dudx[1][j]) < 1e-15)
			uright[0][j] = u[1][j];

		delminus = 2*dudx[N][j]*(x[N+1]-x[N]);
		sminus = lim->limiter_function(delminus, u[N+1][j]-u[N][j]);
		uleft[N][j] = u[N][j] + sminus/4.0*( (1-k*sminus)*delminus + (1+k*sminus)*(u[N+1][j]-u[N][j]) );
			
		if(fabs(dudx[N][j]) < 1e-15)
			uleft[N][j] = u[N][j];
	}
}
