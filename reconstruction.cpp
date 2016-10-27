/** \file reconstruction.cpp
 * \brief Implementation of reconstruction schemes.
 * \author Aditya Kashi
 * \date
 */

#include "reconstruction.hpp"

SlopeReconstruction::SlopeReconstruction(const int _N, const std::vector<double>& _x, const std::vector<double>& _dx, const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: N(_N), x(_x), dx(_dx), u(_u), dudx(_dudx)
{ }
	
SlopeReconstruction::~SlopeReconstruction()
{ }

TrivialSlopeReconstruction::TrivialSlopeReconstruction(const int _N, const std::vector<double>& _x, const std::vector<double>& _dx, 
		const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: SlopeReconstruction(_N, _x, _dx, _u, _dudx)
{ }

void TrivialSlopeReconstruction::compute_slopes()
{
	for(int i = 1; i <= N; i++)
		for(int j = 0; j < NVARS; j++)
			dudx[i][j] = 0;
}

LeastSquaresReconstruction::LeastSquaresReconstruction(const int _N, const std::vector<double>& _x, const std::vector<double>& _dx, 
		const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: SlopeReconstruction(_N, _x, _dx, _u, _dudx)
{ }

void LeastSquaresReconstruction::compute_slopes()
{
	// implement least-squares reconstruction
	double num, denom;
	for(int i = 1; i <= N; i++)
		for(int j = 0; j < NVARS; j++)
		{
			num = (u[i-1][j]-u[i][j])*(x[i-1]-x[i]) + (u[i+1][j]-u[i][j])*(x[i+1]-x[i]);
			denom = (x[i-1]-x[i])*(x[i-1]-x[i]) + (x[i+1]-x[i])*(x[i+1]-x[i]);
			dudx[i][j] = num/denom;
		}

	// one-sided derivatives for ghost cells
	for(int j = 0; j < NVARS; j++)
	{
		dudx[0][j] = (u[1][j] - u[0][j])/(x[1]-x[0]);
		dudx[N+1][j] = (u[N+1][j] - u[N][j])/(x[N+1]-x[N]);
		//dudx[0][j] = 0;
		//dudx[N+1][j] = 0;
	}
}

CentralDifferenceReconstruction::CentralDifferenceReconstruction(const int _N, const std::vector<double>& _x, const std::vector<double>& _dx, 
		const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx)
	: SlopeReconstruction(_N, _x, _dx, _u, _dudx)
{ }

void CentralDifferenceReconstruction::compute_slopes()
{
	for(int i = 1; i <= N; i++)
		for(int j = 0; j < NVARS; j++)
		{
			dudx[i][j] = (u[i+1][j]-u[i-1][j])/(2*dx[i]);
		}

	// one-sided derivatives for ghost cells
	for(int j = 0; j < NVARS; j++)
	{
		dudx[0][j] = (u[1][j] - u[0][j])/(x[1]-x[0]);
		dudx[N+1][j] = (u[N+1][j] - u[N][j])/(x[N+1]-x[N]);
	}
}

TVDSlopeReconstruction::TVDSlopeReconstruction(const int _N, const std::vector<double>& _x, const std::vector<double>& _dx, 
		const std::vector<std::vector<double>>& _u, std::vector<std::vector<double>>& _dudx, std::string _lim)
	: SlopeReconstruction(_N, _x, _dx, _u, _dudx)
{
	if(_lim == "minmod")
	{
		lim = new MinmodLimiter1(SMALL_NUMBER);
		std::cout << "TVDSlopeReconstruction: Using minmod limiter" << std::endl;
	}
	else if(_lim == "vanalbada")
	{
		lim = new VanAlbadaLimiter1(SMALL_NUMBER);
		std::cout << "TVDSlopeReconstruction: Using Van Albada limiter" << std::endl;
	}
	else
	{
		lim = new NoLimiter1(SMALL_NUMBER);
		std::cout << "TVDSlopeReconstruction: Using no limiter!" << std::endl;
	}
}

TVDSlopeReconstruction::~TVDSlopeReconstruction()
{
	delete lim;
}

void TVDSlopeReconstruction::compute_slopes()
{
	int i, j;
	double s1, s2;
	for(i = 1; i <= N; i++)
	{
		for(j = 0; j < NVARS; j++)
		{
			s1 = (u[i+1][j] - u[i][j])/(x[i+1]-x[i]);
			s2 = (u[i][j] - u[i-1][j])/(x[i]-x[i-1]);
			dudx[i][j] = lim->limiter_function(s1,s2);
		}
	}

	// unlimited derivatives for ghost cells
	for(j = 0; j < NVARS; j++)
	{
		dudx[0][j] = (u[1][j]-u[0][j])/(x[1]-x[0]);
		dudx[N+1][j] = (u[N+1][j]-u[N][j])/(x[N+1]-x[N]);
	}
}
	
FaceReconstruction::FaceReconstruction(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, 
		std::vector<std::vector<double>>& _uleft, std::vector<std::vector<double>>& _uright)
	: N(_N), x(_x), u(_u), dudx(_dudx), uleft(_uleft), uright(_uright)
{ }
	
FaceReconstruction::~FaceReconstruction()
{ }

MUSCLReconstruction::MUSCLReconstruction(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, 
		std::vector<std::vector<double>>& uleft, std::vector<std::vector<double>>& uright, std::string _limiter, double _k)
	: FaceReconstruction(_N, _x, _u, _dudx, uleft, uright), limiter(_limiter)
{
	if(limiter == "vanalbada")
	{
		lim = new VanAlbadaLimiter();
		std::cout << "MUSCLReconstruction: Using Van Albada limiter" << std::endl;
	}
	else if(limiter == "minmod")
	{
		lim = new MinmodLimiter();
		std::cout << "MUSCLReconstruction: Using minmod limiter" << std::endl;
	}
	else if(limiter == "hemkerkoren")
	{
		lim = new HemkerKorenLimiter();
		std::cout << "MUSCLReconstruction: Using Hemker-Koren limiter" << std::endl;
	}
	else if(limiter == "vanleer")
	{
		lim = new VanLeerLimiter();
		std::cout << "MUSCLReconstruction: Using Van Leer limiter" << std::endl;
	}
	else
	{
		lim = new NoLimiter();
		std::cout << "MUSCLReconstruction: Caution: not using any limiter.\n";
	}
}

MUSCLReconstruction::~MUSCLReconstruction()
{
	delete lim;
}

void MUSCLReconstruction::compute_face_values()
{
	int i, j, k;
	double denL, denR, num, rL, rR;

	// interior faces
	for(i = 1; i <= N-1; i++)
	{
		for(j = 0; j < NVARS; j++)
		{
			denL = u[i][j]-u[i-1][j];
			denR = u[i+2][j]-u[i+1][j];
			num = u[i+1][j]-u[i][j];

			if(fabs(denL) > ZERO_TOL*10)
			{
				rL = num/denL;
				uleft[i][j] = u[i][j] + 0.25*lim->limiter_function(rL) * ((1.0+k)*num + (1.0-k)*denL);
			}
			else
				uleft[i][j] = u[i][j];

			if(fabs(denR) > ZERO_TOL*10)
			{
				rR = num/denR;
				uright[i][j] = u[i+1][j] - 0.25*lim->limiter_function(rR) * ((1.0+k)*num + (1.0-k)*denR);
			}
			else
				uright[i][j] = u[i+1][j];
		}
	}
	
	// boundaries
	for(j = 0; j < NVARS; j++)
	{
		// extrapolate variables
		double slope, cc, um0, xm0, umN, xmN;

		slope = (u[1][j] - u[0][j])/(x[1]-x[0]);
		cc = u[0][j] - (u[1][j]-u[0][j])/(x[1]-x[0])*x[0];
		xm0 = x[0] - (x[1]-x[0]);
		um0 = slope*xm0 + cc;

		slope = (u[N+1][j] - u[N][j])/(x[N+1]-x[N]);
		cc = u[N][j] - (u[N+1][j]-u[N][j])/(x[N+1]-x[N])*x[N];
		xmN = x[N+1] + x[N]-x[N-1];
		umN = slope*xmN + cc;

		// left
		denL = u[0][j] - um0;
		denR = u[2][j]-u[1][j];
		num = u[1][j]-u[0][j];

		if(fabs(denL) > ZERO_TOL*10)
		{
			rL = num/denL;
			uleft[0][j] = u[0][j] + 0.25*lim->limiter_function(rL) * ((1.0+k)*num + (1.0-k)*denL);
		}
		else
			uleft[0][j] = u[0][j];

		if(fabs(denR) > ZERO_TOL*10)
		{
			rR = num/denR;
			uright[0][j] = u[1][j] - 0.25*lim->limiter_function(rR) * ((1.0+k)*num + (1.0-k)*denR);
		}
		else
			uright[0][j] = u[1][j];
		
		// right
		denL = u[N][j]-u[N-1][j];
		denR = umN - u[N+1][j];
		num = u[N+1][j]-u[N][j];

		if(fabs(denL) > ZERO_TOL*10)
		{
			rL = num/denL;
			uleft[N][j] = u[N][j] + 0.25*lim->limiter_function(rL) * ((1.0+k)*num + (1.0-k)*denL);
		}
		else
			uleft[N][j] = u[N][j];

		if(fabs(denR) > ZERO_TOL*10)
		{
			rR = num/denR;
			uright[N][j] = u[N+1][j] - 0.25*lim->limiter_function(rR) * ((1.0+k)*num + (1.0-k)*denR);
		}
		else
			uright[N][j] = u[N+1][j];
	}
}

MUSCLReconstructionG::MUSCLReconstructionG(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, 
		std::vector<std::vector<double>>& uleft, std::vector<std::vector<double>>& uright, std::string _limiter, double _k)
	: FaceReconstruction(_N, _x, _u, _dudx, uleft, uright), limiter(_limiter)
{
	if(limiter == "vanalbada")
	{
		lim = new VanAlbadaLimiter1(SMALL_NUMBER);
		std::cout << "MUSCLReconstruction: Using Van Albada limiter" << std::endl;
	}
	else if(limiter == "minmod")
	{
		lim = new MinmodLimiter1(SMALL_NUMBER);
		std::cout << "MUSCLReconstruction: Using minmod limiter" << std::endl;
	}
	/*else if(limiter == "hemkerkoren")
	{
		lim = new HemkerKorenLimiter();
		std::cout << "MUSCLReconstruction: Using Hemker-Koren limiter" << std::endl;
	}*/
	else
	{
		lim = new NoLimiter1(SMALL_NUMBER);
		std::cout << "MUSCLReconstruction: Caution: not using any limiter.\n";
	}
}

MUSCLReconstructionG::~MUSCLReconstructionG()
{
	delete lim;
}

void MUSCLReconstructionG::compute_face_values()
{
	/// NOTE: iterates over ALL faces; ghost cells must have valid derivatives!
	int i, j, k;
	double sminus, splus, delminus, delplus;
	for(i = 0; i <= N; i++)
	{
		// get deltas
		for(j = 0; j < NVARS; j++)
		{
			delminus = 2.0*dudx[i][j]*(x[i+1]-x[i]) - (u[i+1][j]-u[i][j]);
			delplus = 2.0*dudx[i+1][j]*(x[i+1]-x[i]) - (u[i+1][j]-u[i][j]);
			
			sminus = lim->limiter_function( delminus, u[i+1][j]-u[i][j] );
			splus = lim->limiter_function( delplus, u[i+1][j]-u[i][j] );

			uleft[i][j] = u[i][j] + sminus/4.0*( (1-k)*delminus + (1+k)*(u[i+1][j]-u[i][j]) );
			uright[i][j] = u[i+1][j] - splus/4.0*( (1-k)*delplus + (1+k)*(u[i+1][j]-u[i][j]) );

			if(fabs(dudx[i][j]) < 1e-15)
				uleft[i][j] = u[i][j];
			if(fabs(dudx[i+1][j]) < 1e-15)
				uright[i][j] = u[i+1][j];

			// hard limit
			if(uleft[i][j] <= ZERO_TOL*100)
				uleft[i][j] = SMALL_NUMBER*100;
			if(uright[i][j] <= ZERO_TOL*100)
				uright[i][j] = SMALL_NUMBER*100;
		}
	}

	// get uright at 0 and uleft at N (boundary faces)
	/*for(j = 0; j < NVARS; j++)
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
	}*/
}

LinearReconstruction::LinearReconstruction(const int _N, const std::vector<double>& _x, const std::vector<std::vector<double>>& _u, const std::vector<std::vector<double>>& _dudx, 
		std::vector<std::vector<double>>& uleft, std::vector<std::vector<double>>& uright)
	: FaceReconstruction(_N, _x, _u, _dudx, uleft, uright)
{
}

void LinearReconstruction::compute_face_values()
{
	int i,j;
	for(i = 0; i <= N; i++)
	{
		for(j = 0; j < NVARS; j++)
		{
			uleft[i][j]  = u[i][j]  +  dudx[i][j]  *  (x[i]-x[i-1])/2.0;
			uright[i][j] = u[i+1][j] + dudx[i+1][j] * (x[i]-x[i+1])/2.0;
		}
	}
}

