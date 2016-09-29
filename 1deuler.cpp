#include "1deuler.hpp"

Euler1d::Euler1d(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs, double CFL, 
		std::string inviscid_flux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter)
	: N(num_cells), domlen(length), bcL(leftBCflag), bcR(rightBCflag), bcvalL(leftBVs), bcvalR(rightBVs), cfl(CFL)
{
	x.resize(N+2);
	dx.resize(N+2);
	u.resize(N+2);
	uleft.resize(N+1);
	uright.resize(N+1);
	dudx.resize(N+2);
	res.resize(N+2);
	A.resize(N+2);
	Af.resize(N+1);
	vol.resize(N+2);
	nodes.resize(N+1);

	for(int i = 0; i < N+2; i++)
	{
		u[i].resize(NVARS);
		dudx[i].resize(NVARS);
		res[i].resize(NVARS);
	}
	for(int i = 0; i < N+1; i++)
	{
		uleft[i].resize(NVARS);
		uright[i].resize(NVARS);
	}

	if(inviscid_flux == "vanleer")
	{
		flux = new VanLeerFlux();
		std::cout << "Euler1d: Using Van Leer numerical flux.\n";
	}
	else if(inviscid_flux == "llf")
	{
		flux = new LocalLaxFriedrichsFlux();
		std::cout << "Euler1d: Using local Lax-Friedrichs numerical flux.\n";
	}

	if(slope_scheme == "none")
	{
		cslope = new TrivialSlopeReconstruction(N,x,u,dudx);
		std::cout << "Euler1d: No slope reconstruction to be used.\n";
	}
	else
	{
		cslope = new LeastSquaresReconstruction(N,x,u,dudx);
		std::cout << "Euler1d: Least-squares slope reconstruction will be used.\n";
	}

	rec = new MUSCLReconstruction(N,x,u,dudx,uleft,uright,limiter,muscl_k);
}

Euler1d::~Euler1d()
{
	delete flux;
	delete cslope;
	delete rec;
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
			A[i] = cellCenteredAreas[i-1];
		/*double h = 0.15, t1 = 0.8, t2 = 3.0;
		for(int i = 1; i < N+1; i++)
			A[i] = 1.0 - h*pow(sin(PI*pow(x[i],t1)),t2);*/

		// maybe assign ghost cell areas by linear extrapolation?
		A[0] = A[1];
		A[N+1] = A[N];
	}

	/** Get interface areas as inverse-distance weighted averages of cell-centered areas so that they are exact for linear profiles.
	 */
	for(int i = 0; i <= N; i++)
		Af[i] = (A[i]*dx[i+1] + A[i+1]*dx[i])/(dx[i]+dx[i+1]);

	for(int i = 0; i < N+2; i++)
	{
		vol[i] = dx[i]*A[i];
	}
}

void Euler1d::compute_inviscid_fluxes()
{
	std::vector<std::vector<double>> fluxes(N+1);

	for(int i = 0; i < N+1; i++)
		fluxes[i].resize(NVARS);

	// iterate over interfaces
	for(int i = 0; i < N+1; i++)
	{
		flux->compute_flux(uleft[i], uright[i], fluxes[i]);

		// update residual
		for(int j = 0; j < NVARS; j++)
		{
			fluxes[i][j] *= Af[i];
			res[i][j] -= fluxes[i][j];
			res[i+1][j] += fluxes[i][j];
		}
	}
	//std::cout << "flux n " << fluxes[N][0] << " " << fluxes[N][1] << std::endl;
}

void Euler1d::compute_source_term()
{
	double p;
	for(int i = 1; i <= N; i++)
	{
		p = (g-1.0)*(u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]);
		res[i][1] += p*(Af[i] - Af[i-1]);
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
		double M_in, c_in, v_in, p_in;
		v_in = u[0][1]/u[0][0];
		p_in = (g-1.0)*(u[0][2] - 0.5*u[0][0]*v_in*v_in);
		c_in = sqrt(g*p_in/u[0][0]);
		M_in = v_in/c_in;

		if(M_in >= 1.0)
		{
			// supersonic inflow
			// get conserved variables from pt, Tt and M
			//double astar = 2*g*(g-1.0)/(g+1.0)*Cv*bcvalL[1];
			double T = bcvalL[1]/(1 + (g-1.0)/2.0*bcvalL[2]*bcvalL[2]);
			double c = sqrt(g*R*T);
			double v = bcvalL[2]*c;
			double p = bcvalL[0]*pow( 1+(g-1.0)/2.0*bcvalL[2]*bcvalL[2], -g/(g-1.0) );
			double rho = p/(R*T);
			double E = p/(g-1.0) + 0.5*rho*v*v;
			/*u[0][0] = 2*rho - u[1][0];
			u[0][1] = 2*rho*v - u[1][1];
			u[0][2] = 2*E - u[1][2];*/
			u[0][0] = rho;
			u[0][1] = rho*v;
			u[0][2] = E;
		}
		else if(M_in >= 0)
		{
			// subsonic inflow
			// get conserved variables from pt and Tt specified in bcvalL[0] and bcvalL[1] respectively
			double vold, pold, cold, vold1, pold1, cold1, astar, dpdu, dt0, lambda, du, v, T, p, c, M;
			std::vector<double> uold0 = u[0];
			vold = uold0[1]/uold0[0];
			pold = (g-1)*(uold0[2] - 0.5*uold0[1]*uold0[1]/uold0[0]);
			cold = sqrt( g*pold/uold0[0] );
			vold1 = u[1][1]/u[1][0];
			pold1 = (g-1)*(u[1][2] - 0.5*u[1][1]*u[1][1]/u[1][0]);
			cold1 = sqrt( g*pold1/u[1][0] );

			astar = 2*g*(g-1.0)/(g+1.0)*Cv*bcvalL[1];
			dpdu = bcvalL[0]*g/(g-1.0)*pow(1.0-(g-1)/(g+1.0)*vold*vold/astar, 1.0/(g-1.0)) * (-2.0)*(g-1)/(g+1.0)*vold/astar;
			dt0 = cfl*dx[0]/(fabs(vold)+cold);
			lambda = (vold1+vold - cold1-cold)*0.5*dt0/dx[0];
			du = -lambda * (pold1-pold-uold0[0]*cold*(vold1-vold)) / (dpdu-uold0[0]*cold);

			v = vold + du;
			T = bcvalL[1]*(1.0 - (g-1)/(g+1.0)*vold*vold/astar);
			p = bcvalL[0]*pow(T/bcvalL[1], g/(g-1.0));
			u[0][0] = p/(R*T);
			u[0][1] = u[0][0]*v;
			u[0][2] = p/(g-1.0) + 0.5*u[0][0]*v*v;

			/*c = sqrt(g*p/u[0][0]);
			M = v/c;
			std::cout << "  apply_boundary_conditions(): Inlet ghost cell mach number = " << M << std::endl;*/
		}
		else
			std::cout << "! Euler1d: apply_boundary_conditions(): Error! Inlet is becoming outlet!" << std::endl;
	}
	else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;

	if(bcR == 0)
	{
		u[N+1][0] = u[N][0];
		u[N+1][1] = -u[N][1];
		u[N+1][2] = u[N][2];
	}
	else if(bcR == 3)
	{
		// outflow
		double l1, l2, l3, cold, cold1, pold, pold1, vold, vold1, dt0, r1, r2, r3, Mold, dp, drho, dv, p;
		std::vector<double> uold = u[N+1];
		vold = uold[1]/uold[0];
		pold = (g-1.0)*(uold[2]-0.5*uold[0]*vold*vold);
		cold = sqrt(g*pold/uold[0]);
		vold1 = u[N][1]/u[N][0];
		pold1 = (g-1.0)*(u[N][2]-0.5*u[N][0]*vold1*vold1);
		cold1 = sqrt(g*pold1/u[N][0]);

		dt0 = cfl*dx[N+1]/(fabs(vold)+cold);
		l1 = (vold+vold1)*0.5*dt0/dx[N+1];
		l2 = (vold+vold1 + cold+cold1)*0.5*dt0/dx[N+1];
		l3 = (vold+vold1 - cold-cold1)*0.5*dt0/dx[N+1];

		r1 = -l1*( u[N+1][0] - u[N][0] - 1.0/(cold*cold)*(pold - pold1));
		r2 = -l2*( pold - pold1 + u[N+1][0]*cold*(vold - vold1));
		r3 = -l3*( pold - pold1 - u[N+1][0]*cold*(vold - vold1));
		Mold = (vold+vold1)/(cold+cold1);

		// check whether supersonic or subsonic
		if(Mold > 1)
			dp = 0.5*(r2+r3);
		else
			dp = 0;

		drho = r1 + dp/(cold*cold);
		dv = (r2-dp)/(u[N+1][0]*cold);

		u[N+1][0] += drho;
		u[N+1][1] = u[N+1][0]*(vold + dv);

		if(Mold > 1)
			p = pold + dp;
		else
			p = bcvalR[0];

		u[N+1][2] = p/(g-1.0) + 0.5*u[N+1][0]*(vold+dv)*(vold+dv);
	}
	else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;
}
	
void Euler1d::apply_boundary_conditions_at_left_boundary(std::vector<double>& ul, const std::vector<double>& ur)
{
	if(bcL == 0)
	{
		// homogeneous Dirichlet - wall
		ul[0] = ur[0];
		ul[1] = -ur[1];
		ul[2] = ur[2];
	}
	else if(bcL == 1)
	{
		double M_in, c_in, v_in, p_in;
		v_in = ur[1]/ur[0];
		p_in = (g-1.0)*(ur[2] - 0.5*ur[0]*v_in*v_in);
		c_in = sqrt(g*p_in/ur[0]);
		M_in = v_in/c_in;
		
		/*double M_ex, c_ex, v_ex, p_ex;
		v_ex = ul[1]/ul[0];
		p_ex = (g-1.0)*(ul[2] - 0.5*ul[0]*v_ex*v_ex);
		c_ex = sqrt(g*p_ex/ul[0]);
		M_ex = v_ex/c_ex;

		double M_eff = (M_in+M_ex)/2.0;*/
		
		// get fluid state from presribed free-stream conditions
		double T = bcvalL[1]/(1 + (g-1.0)/2.0*bcvalL[2]*bcvalL[2]);
		double c = sqrt(g*R*T);
		double v = bcvalL[2]*c;
		double p = bcvalL[0]*pow( 1+(g-1.0)/2.0*bcvalL[2]*bcvalL[2], -g/(g-1.0) );
		double rho = p/(R*T);
		double E = p/(g-1.0) + 0.5*rho*v*v;

		double M_eff = (M_in+bcvalL[2])*0.5;

		if(M_eff >= 1.0)
		{
			/// supersonic inflow
			/// get conserved variables from prescribed pt, Tt and M
			ul[0] = rho;
			ul[1] = rho*v;
			ul[2] = E;
		}
		else if(M_eff >= 0)
		{
			/// Subsonic inflow
			/** Get Riemann invarients corresponding to eigenvalues u-c, u and u+c resp.
			 * Take R1 from interior, and R2 and R3 from free stream
			 */
			double R1, R2, R3, vg, cg, pg;
			R1 = v_in - 2.0*c_in/(g-1.0);
			R2 = p/pow(rho,g);
			R3 = v + 2.0*c/(g-1.0);
			vg = (R3+R1)*0.5;
			cg = (R3-R1)*(g-1.0)/4.0;
			ul[0] = pow(cg*cg/(g*R2), 1.0/(g-1.0));
			pg = ul[0]*cg*cg/g;
			ul[1] = ul[0]*vg;
			ul[2] = pg/(g-1.0) + 0.5*ul[0]*vg*vg;
		}
		else
			std::cout << "! Euler1d: apply_boundary_conditions_left(): Error! Inlet is becoming outlet!" << std::endl;
	}
	else std::cout << "! Euler1D: apply_boundary_conditions_left(): BC type not recognized!" << std::endl;
}

void Euler1d::apply_boundary_conditions_at_right_boundary(const std::vector<double>& ul, std::vector<double>& ur)
{
	if(bcR == 0)
	{
		ur[0] = ul[0];
		ur[1] = -ul[1];
		ur[2] = ul[2];
	}
	else if(bcR == 3)
	{
		// outflow
		double R1, R2, R3, cold0, pold0, vold0, cold1, pold1, vold, vold1, Minf, pinf, Tinf, rhoinf, vinf, cinf, v, c, p, Meff;
		pinf = bcvalR[0]; Tinf = bcvalR[1]; Minf = bcvalR[2];
		rhoinf = pinf/(R*Tinf);
		cinf = sqrt(g*pinf/rhoinf);
		vinf = Minf*cinf;

		vold1 = ul[1]/ul[0];
		pold1 = (g-1.0)*(ul[2]-0.5*ul[0]*vold1*vold1);
		cold1 = sqrt(g*pold1/ul[0]);

		// previous values at the ghost point we want to set
		vold0 = ur[1]/ur[0];
		pold0 = (g-1.0)*(ur[2]-0.5*ur[0]*vold0*vold0);
		cold0 = sqrt(g*pold0/ur[0]);

		Meff = 0.5*(vold1/cold1 + Minf);
		
		if(Meff > 0 && Meff < 1.0)
			R1 = vinf - 2*cinf/(g-1.0);
		else if(Meff >= 1)
			R1 = vold1 - 2*cold1/(g-1.0);
		else
			std::cout << "! Euler1d: apply_boundary_conditions_at_right_boundary(): Flow is negative!" << std::endl;

		R2 = pold1/pow(ul[0],g);
		R3 = vold1 + 2*cold1/(g-1.0);

		v = 0.5*(R3+R1); c = 0.25*(g-1.0)*(R3-R1);
		ur[0] = pow(c*c/(g*R2), 1.0/(g-1.0));
		p = ur[0]*c*c/g;
		ur[1] = ur[0]*v;
		ur[2] = p/(g-1.0) + 0.5*ur[0]*v*v;
		
		/*R3 = vold1 + 2*cold1/(g-1.0);

		v = 0.5*(R3+R1); c = 0.25*(g-1.0)*(R3-R1);
		p = pinf;
		ur[0] = g*p/(c*c);
		ur[1] = ur[0]*v;
		ur[2] = p/(g-1.0) + 0.5*ur[0]*v*v;*/
	}
	else std::cout << "! Euler1D: apply_boundary_conditions_right(): BC type not recognized!" << std::endl;
}


Euler1dExplicit::Euler1dExplicit(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs,
		double CFL, std::string inviscid_flux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter, double f_time, int temporal_order, std::string RKfile)
	: Euler1d(num_cells,length,leftBCflag,rightBCflag,leftBVs,rightBVs, CFL, inviscid_flux, slope_scheme, face_extrap_scheme, limiter), ftime(f_time), temporalOrder(temporal_order)
{
	maxWaveSpeed.resize(N+2);
	RKCoeffs.resize(temporalOrder);
	for(int i = 0; i < temporalOrder; i++)
		RKCoeffs[i].resize(3);

	std::ifstream rkfile(RKfile);
	for(int i = 0; i < temporalOrder; i++)
		for(int j = 0; j < 3; j++)
			rkfile >> RKCoeffs[i][j];
	rkfile.close();

	std::cout << "Euler1dExplicit: Using " << temporalOrder << "-stage TVD RK scheme; loaded coefficients.\n";
}

void Euler1dExplicit::run()
{
	int step = 0, istage;
	double dt, time = 0;

	// IC for Sod shock tube
	for(int i = 0; i < N+2; i++)
		if(x[i] <= 0.5)
		{
			u[i][0] = 1.0;
			u[i][1] = 0;
			u[i][2] = 2.5;
		}
		else
		{
			u[i][0] = 0.125;
			u[i][1] = 0;
			u[i][2] = 0.25;
		}

	std::vector<double> c(N+2);
	std::vector<std::vector<double>> uold, ustage;
	uold.resize(N+2);
	ustage.resize(N+2);
	for(int i = 0; i < N+2; i++)
	{
		uold[i].resize(NVARS);
		ustage[i].resize(NVARS);
	}

	while(time < ftime)
	{
		int i,j;
		for(i = 0; i < N+2; i++)
		{
			for(j = 0; j < NVARS; j++)
			{
				uold[i][j] = u[i][j];
				res[i][j] = 0;
			}
		}
		
		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		for(i = 1; i < N+1; i++)
		{
			c[i] = sqrt( g*(g-1.0) * (u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]) / u[i][0] );
		}

		double mws = dx[1]/(fabs(u[1][1]) + c[1]);
		double a;
		for(i = 2; i < N+1; i++)
		{
			a = dx[i]/(fabs(u[i][1]) + c[i]);
			if(a < mws) mws = a;
		}

		dt = cfl*mws;

		for(istage = 0; istage < temporalOrder; istage++)
		{
			for(i = 0; i < N+2; i++)
			{
				for(j = 0; j < NVARS; j++)
				{
					ustage[i][j] = u[i][j];
					res[i][j] = 0;
				}
			}
			cslope->compute_slopes();
			rec->compute_face_values();
			/*apply_boundary_conditions_at_left_boundary(uleft[0], uright[0]);
			apply_boundary_conditions_at_right_boundary(uleft[N], uright[N]);*/

			compute_inviscid_fluxes();
			compute_source_term();

			// RK stage
			for(i = 1; i < N+1; i++)
				for(j = 0; j < NVARS; j++)
					u[i][j] = RKCoeffs[istage][0]*uold[i][j] + RKCoeffs[istage][1]*ustage[i][j] + RKCoeffs[istage][2]*dt/vol[i]*res[i][j];

			// apply BCs
			/*apply_boundary_conditions_at_left_boundary(u[0], u[1]);
			apply_boundary_conditions_at_right_boundary(u[N], u[N+1]);*/
			apply_boundary_conditions();
		}

		if(step % 10 == 0)
			std::cout << "Euler1dExplicit: run(): Step " << step << " - Time = " << time << std::endl;

		time += dt;
		step++;
	}

	std::cout << "Euler1dExplicit: run(): Done. Number of time steps = " << step << ", final time = " << time << std::endl;
}

void Euler1dExplicit::postprocess(std::string outfilename)
{
	std::ofstream ofile(outfilename);
	double pressure, mach, c;
	for(int i = 1; i < N+1; i++)
	{
		pressure = (g-1)*(u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]);
		c = sqrt(g*pressure/u[i][0]);
		mach = (u[i][1]/u[i][0])/c;
		ofile << x[i] << " " << u[i][0] << " " << mach << " " << pressure << " " << u[i][1] << " " << c << '\n';
	}
	ofile.close();
}


Euler1dSteadyExplicit::Euler1dSteadyExplicit(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs, 
		double CFL, std::string inviscidFlux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter, double toler, int max_iter)
	: Euler1d(num_cells,length,leftBCflag,rightBCflag,leftBVs,rightBVs, CFL, inviscidFlux, slope_scheme, face_extrap_scheme, limiter), tol(toler), maxiter(max_iter)
{
	maxWaveSpeed.resize(N+2);
}

void Euler1dSteadyExplicit::run()
{
	int step = 0;
	double resnorm = 1.0, resnorm0 = 1.0;
	std::vector<double> dt(N+2);

	std::vector<double> c(N+2);
	std::vector<std::vector<double>> uold;
	uold.resize(N+2);
	for(int i = 0; i < N+2; i++)
		uold[i].resize(NVARS);

	// initial conditions
	
	double p_t = bcvalL[0];
	double T_t = bcvalL[1];
	double M = bcvalL[2];
	double term = 1.0 + (g-1.0)*0.5*M*M;
	double Tin = T_t/term;
	double pin = p_t*pow(term, -g/(g-1.0));

	double pex = bcvalR[0], tex = bcvalR[1], Mex = bcvalR[2], cex, vex;

	double cin;
	
	// set some cells according to inlet condition
	for(int i = 0; i <= N; i++)
	{
		u[i][0] = pin/(R*Tin);
		double cin = sqrt(g*pin/u[i][0]);
		u[i][1] = u[i][0]*M*cin;
		u[i][2] = pin/(g-1.0)+0.5*u[i][1]*M*cin;
	}

	// set rest of the cells according to exit conditions
	for(int i = N+1; i <= N+1; i++)
	{
		u[i][0] = pex/(R*tex);
		cex = sqrt(g*pex/u[i][0]);
		vex = Mex*cex;
		u[i][1] = u[i][0]*vex;
		u[i][2] = pex/(g-1.0) + 0.5*u[i][0]*vex*vex;
	}

	// Start time loop

	while(resnorm/resnorm0 > tol && step < maxiter)
	{
		int i,j;
		for(i = 0; i < N+2; i++)
		{
			for(j = 0; j < NVARS; j++)
			{
				uold[i][j] = u[i][j];
				res[i][j] = 0;
			}
		}

		cslope->compute_slopes();
		rec->compute_face_values();
		/*apply_boundary_conditions_at_left_boundary(uleft[0], uright[0]);
		apply_boundary_conditions_at_right_boundary(uleft[N], uright[N]);*/

		compute_inviscid_fluxes();
		compute_source_term();

		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		for(i = 1; i < N+1; i++)
		{
			c[i] = sqrt( g*(g-1.0) * (u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]) / u[i][0] );
		}

		for(i = 1; i < N+1; i++)
		{
			dt[i] = cfl * dx[i]/(fabs(u[i][1]) + c[i]);
		}

		resnorm = 0;
		for(i = 1; i <= N; i++)
			resnorm += res[i][0]*res[i][0]*dx[i];
		resnorm = sqrt(resnorm);
		if(step==0)
			resnorm0 = resnorm;

		// RK step
		for(i = 1; i < N+1; i++)
			for(j = 0; j < NVARS; j++)
				u[i][j] = uold[i][j] + dt[i]/vol[i]*res[i][j];

		// apply BCs
		/*apply_boundary_conditions_at_left_boundary(u[0], u[1]);
		apply_boundary_conditions_at_right_boundary(u[N], u[N+1]);*/
		apply_boundary_conditions();

		if(step % 10 == 0)
			std::cout << "Euler1dSteadyExplicit: run(): Step " << step << ", relative mass flux norm = " << resnorm/resnorm0 << std::endl;

		step++;
	}

	std::cout << "Euler1dExplicit: run(): Done. Number of time steps = " << step << std::endl;
	if(step == maxiter)
		std::cout << "Euler1dExplicit: run(): Not converged!" << std::endl;
}

void Euler1dSteadyExplicit::postprocess(std::string outfilename)
{
	std::ofstream ofile(outfilename);
	double pressure, mach, c;
	for(int i = 1; i < N+1; i++)
	{
		pressure = (g-1)*(u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]);
		c = sqrt(g*pressure/u[i][0]);
		mach = (u[i][1]/u[i][0])/c;
		ofile << x[i] << " " << u[i][0] << " " << mach << " " << pressure/bcvalL[0] << " " << u[i][1] << " " << c << '\n';
	}
	ofile.close();
}
