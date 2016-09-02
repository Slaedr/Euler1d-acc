#include "inviscidflux.hpp"

VanLeerFlux::VanLeerFlux()
{
	fluxL.resize(NVARS);
	fluxR.resize(NVARS);
}

void VanLeerFlux::compute_flux(const std::vector<double>* uleft, const std::vector<double>* uright, std::vector<double>* flux)
{
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

	for(j = 0; j < NVARS; j++)
		flux[j] = fluxL[j] + fluxR[j];
}
