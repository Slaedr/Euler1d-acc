#include "1deuler.hpp"

Euler1d::Euler1d(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs, std::string inviscid_flux)
	: N(num_cells), domlen(length), bcL(leftBCflag), bcR(rightBCflag), bcvalL(leftBVs), bcvalR(rightBVs), inviscidflux(inviscid_flux)
{
	x.resize(N+2);
	dx.resize(N+2);
	u.resize(N+2);
	res.resize(N+2);
	A.resize(N+2);
	vol.resize(N+2);
	nodes.resize(N+1);

	if(inviscidflux == "vanleer")
		flux = new VanLeerFlux();
	else if(inviscidflux == "llf")
		flux = new LocalLaxFriedrichsFlux();
}

Euler1d::~Euler1d()
{
	delete flux;
}

void Euler1d::generate_mesh(int type, const std::vector<double>& pointlist)
{
	if(type == 0)
	{
		nodes[0] = 0.0;
		double delx = domlen/N;
		x[0] = -delx/2.0;
		dx[0] = delx;
		for(int i = 1; i < N+2; i++)
		{
			dx[i] = delx;
			x[i] = i*delx - delx/2.0;
			if(i < N+1)
				nodes[i] = i*delx;
		}
	}
	else
	{
		for(int i = 0; i < N+1; i++)
			nodes[i] = pointlist[i];
		for(int i = 1; i < N+1; i++)
		{
			x[i] = (nodes[i]+nodes[i-1])/2.0;
			dx[i] = nodes[i]-nodes[i-1];
		}

		// ghost cells
		x[0] = -x[1]; 
		dx[0] = dx[1];
		x[N+1] = nodes[N] + dx[N]/2.0;
		dx[N+1] = dx[N];
	}
}

void Euler1d::set_area(int type, std::vector<double>& cellCenteredAreas)
{
	if(type == 0)
		for(int i = 0; i < N+2; i++)
			A[i] = cellCenteredAreas[0];
	else
	{
		for(int i = 1; i < N+1; i++)
			A[i] = cellCenteredAreas[i];

		// TODO: assign ghost cell areas by linear extrapolation
		A[0] = A[1];
		A[N+1] = A[N];
	}

	for(int i = 0; i < N+2; i++)
		vol[i] = dx[i]*A[i];
}

void Euler1d::compute_inviscid_fluxes()
{
	std::vector<std::vector<double>> fluxes(N+1);
	double farea;

	// iterate over interfaces
	for(int i = 0; i < N+1; i++)
	{
		fluxes[i].resize(NVARS);
		flux->compute_flux(u[i], u[i+1], fluxes[i]);

		// get cross-sectional area at the ith section as average of cell-centered values
		farea = (A[i] + A[i+1])/2.0;

		// update residual
		for(int j = 0; j < NVARS; j++)
		{
			fluxes[i][j] *= farea;
			res[i][j] -= vol[i]*fluxes[i][j];
			res[i+1][j] += vol[i+1]*fluxes[i][j];
		}
	}
}

void Euler1d::compute_source_term()
{
	double p;
	for(int i = 1; i <= N; i++)
	{
		p = (g-1.0)*(u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]);
		res[i][1] += p*(A[i+1] - 2*A[i] + A[i-1]);
	}
}

void Euler1d::apply_boundary_conditions()
{
	if(bcL == 0)
	{
		// homogeneous Dirichlet - wall
		u[0][0] = u[1][0];
		u[0][1] = -u[1][1];
		u[0][2] = u[1][2];
	}
	else if(bcL == 1)
	{
		// TODO: supersonic inflow
		u[0][0] = 2*bcvalL[0]-u[1][0];
	}
	else if(bcL == 2)
	{
		// TODO: subsonic inflow
	}
	else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;

	if(bcR == 0)
	{
		u[N+1][0] = u[N][0];
		u[N+1][1] = -u[N+1][1];
		u[N+1][2] = u[N+1][2];
	}
	else if(bcR == 3)
	{
		// TODO: subsonic outflow
	}
	else if(bcR == 4)
	{
		// TODO: supersonic outflow
	}
	else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;
}

Euler1dExplicit::Euler1dExplicit(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs, std::string inviscid_flux, 
		double CFL, double f_time, int temporal_order)
	: Euler1d(num_cells,length,leftBCflag,rightBCflag,leftBVs,rightBVs, inviscid_flux), cfl(CFL), ftime(f_time), temporalOrder(temporalOrder)
{
	maxWaveSpeed.resize(N+2);
}

void Euler1dExplicit::run()
{
	int step;
	double dt, time = 0;

	std::vector<double> c(N+2);
	std::vector<std::vector<double>> uold[2];
	uold[0].resize(N+2);
	for(int i = 0; i < N+2; i++)
		uold[0][i].reserve(NVARS);
	if(temporalOrder > 1)
	{
		uold[1].resize(N+2);
		for(int i = 0; i < N+2; i++)
			uold[1][i].reserve(NVARS);
	}

	while(time < ftime)
	{
		int i,j;
		for(i = 0; i < N+2; i++)
		{
			for(j = 0; j < NVARS; j++)
			{
				uold[0][i][j] = u[i][j];
				res[i][j] = 0;
			}
		}

		compute_inviscid_fluxes();
		compute_source_term();

		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		for(i = 1; i < N+1; i++)
			c[i] = sqrt( g*(g-1.0) * (u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]) / u[i][0] );

		double mws = dx[1]/fabs(u[1][1]) + c[1];
		double a;
		for(i = 2; i < N+1; i++)
		{
			a = dx[i]/fabs(u[i][1]) + c[i];
			if(a < mws) mws = a;
		}

		dt = cfl*mws;

		// RK step
		for(i = 1; i < N+1; i++)
			for(j = 0; j < NVARS; j++)
				u[i][j] = uold[0][i][j] + dt*res[i][j];

		// apply BCs
		apply_boundary_conditions();

		if(step % 10 == 0)
			std::cout << "Euler1dExplicit: run(): Step " << step << " - Time = " << time << std::endl;

		time += dt;
		step++;
	}
}
