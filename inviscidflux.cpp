#include "inviscidflux.hpp"

InviscidFlux::InviscidFlux()
{
	#pragma acc enter data copyin(this)
}

InviscidFlux::~InviscidFlux()
{
	#pragma acc exit data delete(this)
}

LocalLaxFriedrichsFlux::LocalLaxFriedrichsFlux() : InviscidFlux()
{
}

void LocalLaxFriedrichsFlux::compute_flux(double const *const uleft, double const *const uright, double *const flux)
{
	double eps = 0.5;

	double pl = (g-1) * (uleft[2]-0.5*uleft[1]*uleft[1]/uleft[0]);
	double pr = (g-1) * (uright[2]-0.5*uright[1]*uright[1]/uright[0]);
	double cl = sqrt(g* pl / uleft[0]);
	double cr = sqrt(g* pr / uright[0]);

	// find max abs eigenvalue at face
	double el = fabs(uleft[1]/uleft[0]) + cl;
	double er = fabs(uright[1]/uright[0]) + cr;
	double emax = el > er ? el : er;

	// get flux
	flux[0] = 0.5*( uleft[1] + uright[1] - eps*emax*(uright[0]-uleft[0]) );
	flux[1] = 0.5*( uleft[1]*uleft[1]/uleft[0] + pl + uright[1]*uright[1]/uright[0] + pr - eps*emax*(uright[1]-uleft[1]) );
	flux[2] = 0.5*( uleft[1]/uleft[0]*(uleft[2]+pl) + uright[1]/uright[0]*(uright[2]+pr) - eps*emax*(uright[2]-uleft[2]) );
}

void LocalLaxFriedrichsFlux::compute_flux_prim(double const *const uleft, double const *const uright, double *const flux)
{
	double eps = 0.5;

	/*if(uleft[2] < SMALL_NUMBER)
		std::cout << "LLF: " << uleft[2] << std::endl;
	if(uleft[0] < SMALL_NUMBER)
		std::cout << "LLF: " << uleft[0] << std::endl;*/

	double cl = sqrt(g* uleft[2] / uleft[0]);
	double cr = sqrt(g* uright[2] / uright[0]);
	double El = uleft[2]/(g-1.0) + 0.5*uleft[0]*uleft[1]*uleft[1];
	double Er = uright[2]/(g-1.0) + 0.5*uright[0]*uright[1]*uright[1];

	// find max abs eigenvalue at face
	double el = fabs(uleft[1]) + cl;
	double er = fabs(uright[1]) + cr;
	double emax = el > er ? el : er;

	// get flux
	flux[0] = 0.5*( uleft[0]*uleft[1] + uright[0]*uright[1] - eps*emax*(uright[0]-uleft[0]) );
	flux[1] = 0.5*( uleft[1]*uleft[1]*uleft[0] + uleft[2] + uright[1]*uright[1]*uright[0] + uright[2] - eps*emax*(uright[0]*uright[1]-uleft[0]*uleft[1]) );
	flux[2] = 0.5*( uleft[1]*(El+uleft[2]) + uright[1]*(Er+uright[2]) - eps*emax*(Er-El) );
}

VanLeerFlux::VanLeerFlux() : InviscidFlux()
{
}

void VanLeerFlux::compute_flux(double const *const uleft, double const *const uright, double *const flux)
{
	double fluxL[NVARS];
	double fluxR[NVARS];
	int j;
	double Mi, Mj, pi, pj, ci, cj, vi, vj;
	vi = uleft[1]/uleft[0];
	pi = (g-1.0)*(uleft[2] - 0.5*uleft[0]*vi*vi);
	ci = sqrt(g*pi/uleft[0]);
	Mi = vi/ci;
	
	vj = uright[1]/uright[0];
	pj = (g-1.0)*(uright[2] - 0.5*uright[0]*vj*vj);
	cj = sqrt(g*pj/uright[0]);
	Mj = vj/cj;

	if(Mi <= -1)
		for(j = 0; j < NVARS; j++)
			fluxL[j] = 0;
	else if(Mi <= 1)
	{
		fluxL[0] = 0.25*uleft[0]*ci*(Mi+1)*(Mi+1);
		fluxL[1] = fluxL[0]* ((g-1)*Mi + 2)*ci/g;
		fluxL[2] = fluxL[0]* (((g-1)*Mi + 2)*ci)*(((g-1)*Mi + 2)*ci)/(2*(g*g-1.0));
	}
	else if(Mi > 1.0)
	{
		fluxL[0] = uleft[1];
		fluxL[1] = uleft[1]*vi + pi;
		fluxL[2] = vi*(uleft[2] + pi);
	}

	if(Mj >= 1.0)
		for(j = 0; j < NVARS; j++)
			fluxR[j] = 0;
	else if(Mj >= -1.0)
	{
		fluxR[0] = -0.25*uright[0]*cj*(Mj-1)*(Mj-1);
		fluxR[1] = fluxR[0]* ((g-1)*Mj - 2)*cj/g;
		fluxR[2] = fluxR[0]* (((g-1)*Mj - 2)*cj)*(((g-1)*Mj - 2)*cj)/(2*(g*g-1.0));
	}
	else
	{
		fluxR[0] = uright[1];
		fluxR[1] = uright[1]*vj + pj;
		fluxR[2] = vj*(uright[2] + pj);
	}

	for(j = 0; j < NVARS; j++)
		flux[j] = fluxL[j] + fluxR[j];
}

void VanLeerFlux::compute_flux_prim(double const *const uleft, double const *const uright, double *const flux)
{
	double fluxL[NVARS];
	double fluxR[NVARS];
	int j;
	double Mi, Mj, ci, cj;
	ci = sqrt(g*uleft[2]/uleft[0]);
	Mi = uleft[1]/ci;
	
	cj = sqrt(g*uright[2]/uright[0]);
	Mj = uright[1]/cj;

	if(Mi <= -1)
		for(j = 0; j < NVARS; j++)
			fluxL[j] = 0;
	else if(Mi <= 1)
	{
		fluxL[0] = 0.25*uleft[0]*ci*(Mi+1)*(Mi+1);
		fluxL[1] = fluxL[0]* ((g-1)*Mi + 2)*ci/g;
		fluxL[2] = fluxL[0]* (((g-1)*Mi + 2)*ci)*(((g-1)*Mi + 2)*ci)/(2*(g*g-1.0));
	}
	else if(Mi > 1.0)
	{
		fluxL[0] = uleft[0]*uleft[1];
		fluxL[1] = uleft[0]*uleft[1]*uleft[1] + uleft[2];
		fluxL[2] = uleft[1]*(uleft[2]/(g-1.0)+0.5*uleft[0]*uleft[1]*uleft[1] + uleft[2]);
	}

	if(Mj >= 1.0)
		for(j = 0; j < NVARS; j++)
			fluxR[j] = 0;
	else if(Mj >= -1.0)
	{
		fluxR[0] = -0.25*uright[0]*cj*(Mj-1)*(Mj-1);
		fluxR[1] = fluxR[0]* ((g-1)*Mj - 2)*cj/g;
		fluxR[2] = fluxR[0]* (((g-1)*Mj - 2)*cj)*(((g-1)*Mj - 2)*cj)/(2*(g*g-1.0));
	}
	else
	{
		fluxR[0] = uright[0]*uright[1];
		fluxR[1] = uright[0]*uright[1]*uright[1] + uright[2];
		fluxR[2] = uright[1]*(uright[2]/(g-1.0)+0.5*uright[0]*uright[1]*uright[1] + uright[2]);
	}

	for(j = 0; j < NVARS; j++)
		flux[j] = fluxL[j] + fluxR[j];
}

void compute_vanleerflux_prim(double const *const uleft, double const *const uright, double *const flux, const double g)
{
	double fluxL[NVARS];
	double fluxR[NVARS];
	int j;
	double Mi, Mj, ci, cj;
	ci = sqrt(g*uleft[2]/uleft[0]);
	Mi = uleft[1]/ci;
	
	cj = sqrt(g*uright[2]/uright[0]);
	Mj = uright[1]/cj;

	if(Mi <= -1)
		for(j = 0; j < NVARS; j++)
			fluxL[j] = 0;
	else if(Mi <= 1)
	{
		fluxL[0] = 0.25*uleft[0]*ci*(Mi+1)*(Mi+1);
		fluxL[1] = fluxL[0]* ((g-1)*Mi + 2)*ci/g;
		fluxL[2] = fluxL[0]* (((g-1)*Mi + 2)*ci)*(((g-1)*Mi + 2)*ci)/(2*(g*g-1.0));
	}
	else if(Mi > 1.0)
	{
		fluxL[0] = uleft[0]*uleft[1];
		fluxL[1] = uleft[0]*uleft[1]*uleft[1] + uleft[2];
		fluxL[2] = uleft[1]*(uleft[2]/(g-1.0)+0.5*uleft[0]*uleft[1]*uleft[1] + uleft[2]);
	}

	if(Mj >= 1.0)
		for(j = 0; j < NVARS; j++)
			fluxR[j] = 0;
	else if(Mj >= -1.0)
	{
		fluxR[0] = -0.25*uright[0]*cj*(Mj-1)*(Mj-1);
		fluxR[1] = fluxR[0]* ((g-1)*Mj - 2)*cj/g;
		fluxR[2] = fluxR[0]* (((g-1)*Mj - 2)*cj)*(((g-1)*Mj - 2)*cj)/(2*(g*g-1.0));
	}
	else
	{
		fluxR[0] = uright[0]*uright[1];
		fluxR[1] = uright[0]*uright[1]*uright[1] + uright[2];
		fluxR[2] = uright[1]*(uright[2]/(g-1.0)+0.5*uright[0]*uright[1]*uright[1] + uright[2]);
	}

	for(j = 0; j < NVARS; j++)
		flux[j] = fluxL[j] + fluxR[j];
}
