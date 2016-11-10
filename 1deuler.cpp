#include "1deuler.hpp"

Euler1d::Euler1d(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs, double CFL, 
		std::string inviscid_flux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter)
	: N(num_cells), domlen(length), bcL(leftBCflag), bcR(rightBCflag), cfl(CFL)
{
	ncell = N+2;
	nface = N+1;
	x = (double*)malloc((N+2)*sizeof(double));
	dx = (double*)malloc((N+2)*sizeof(double));
	A = (double*)malloc((N+2)*sizeof(double));
	Af = (double*)malloc((N+1)*sizeof(double));
	vol = (double*)malloc((N+2)*sizeof(double));
	nodes = (double*)malloc((N+1)*sizeof(double));

	u = (double**)malloc((N+2)*sizeof(double*));
	//u[0] = (double*)malloc((N+2)*NVARS*sizeof(double));
	prim = (double**)malloc((N+2)*sizeof(double*));
	prim[0] = (double*)malloc((N+2)*NVARS*sizeof(double));
	dudx = (double**)malloc((N+2)*sizeof(double*));
	dudx[0] = (double*)malloc((N+2)*NVARS*sizeof(double));
	res = (double**)malloc((N+2)*sizeof(double*));
	res[0] = (double*)malloc((N+2)*NVARS*sizeof(double));

	/*uleft = (double**)malloc((N+1)*sizeof(double*));
	uleft[0] = (double*)malloc((N+1)*NVARS*sizeof(double));
	uright = (double**)malloc((N+1)*sizeof(double*));
	uright[0] = (double*)malloc((N+1)*NVARS*sizeof(double));*/
	prleft = (double**)malloc((N+1)*sizeof(double*));
	prleft[0] = (double*)malloc((N+1)*NVARS*sizeof(double));
	prright = (double**)malloc((N+1)*sizeof(double*));
	prright[0] = (double*)malloc((N+1)*NVARS*sizeof(double));

	for(int i = 0; i < NVARS; i++)
	{
		bcvalL[i] = leftBVs[i];
		bcvalR[i] = rightBVs[i];
	}

	for(int i = 0; i < ncell; i++)
		u[i] = (double*)malloc(NVARS*sizeof(double));

	for(int i = 1; i < N+2; i++)
	{
		//u[i] = *u + i*NVARS;
		prim[i] = *prim + i*NVARS;
		dudx[i] = *dudx + i*NVARS;
		res[i] = *res + i*NVARS;
	}
	for(int i = 0; i < N+1; i++)
	{
		/*uleft[i] = *uleft + i*NVARS;
		uright[i] = *uright + i*NVARS;*/
		prleft[i] = *prleft + i*NVARS;
		prright[i] = *prright + i*NVARS;
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
		cslope = new TrivialSlopeReconstruction(N,x,dx,u,dudx);
		std::cout << "Euler1d: No slope reconstruction to be used.\n";
	}
	else if(slope_scheme == "leastsquares")
	{
		cslope = new LeastSquaresReconstruction(N,x,dx,prim,dudx);
		std::cout << "Euler1d: Least-squares slope reconstruction will be used.\n";
	}
	else if(slope_scheme == "central")
	{
		cslope = new CentralDifferenceReconstruction(N,x,dx,prim,dudx);
		std::cout << "Euler1d: Central difference slope will be used.\n";
	}
	else
	{
		cslope = new TVDSlopeReconstruction(N,x,dx,prim,dudx,limiter);
		std::cout << "Euler1d: TVD slope reconstruction will be used.\n";
	}

	if(face_extrap_scheme == "MUSCL")
	{
		rec = new MUSCLReconstruction(N,x,prim,dudx,prleft,prright,limiter,muscl_k);
		std::cout << "Euler1d: Using MUSCL face reconstruction." << std::endl;
	}
	else if(face_extrap_scheme == "MUSCLG")
	{
		rec = new MUSCLReconstructionG(N,x,prim,dudx,prleft,prright,limiter,muscl_k);
		std::cout << "Euler1d: Using 'MUSCL-G' face reconstruction." << std::endl;
	}
	else
	{
		rec = new LinearReconstruction(N,x,prim,dudx,prleft,prright);
		std::cout << "Euler1d: Using Linear Taylor expansion reconstruction." << std::endl;
	}

}

Euler1d::~Euler1d()
{
	delete flux;
	delete cslope;
	delete rec;
	
	free(x);			
	free(dx);			
	free(vol);		
	free(nodes);		
	free(A);
	free(Af);

	//free(u[0]);			
	free(prim[0]);		
	//free(uleft[0]);		
	//free(uright[0]);	
	free(prleft[0]);	
	free(prright[0]);	
	free(dudx[0]);		
	free(res[0]);

	for(int i = 0; i < N+2; i++)
		free(u[i]);
	free(u);
	free(prim);		
	//free(uleft);		
	//free(uright);	
	free(prleft);	
	free(prright);	
	free(dudx);		
	free(res);
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

void Euler1d::compute_inviscid_fluxes(double** prleft, double** prright, double** res, double* Af)
{
	double** fluxes = (double**)malloc((N+1)*sizeof(double*));
	fluxes[0] = (double*)malloc(NVARS*(N+1)*sizeof(double));
	for(int i = 0; i < N+1; i++)
		fluxes[i] = *fluxes + i*NVARS;

	#pragma acc kernels present( prleft, prright, Af, res) create(fluxes[:N+1][:NVARS])
	{
		// iterate over interfaces
		#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 0; i < N+1; i++)
		{
			//flux->compute_flux_prim(prleft[i], prright[i], fluxes[i]);
			compute_vanleerflux_prim(prleft[i], prright[i], fluxes[i]);

			// update residual
			for(int j = 0; j < NVARS; j++)
			{
				fluxes[i][j] *= Af[i];
				res[i][j] -= fluxes[i][j];
				res[i+1][j] += fluxes[i][j];
			}
		}
	}
	
	free(fluxes[0]);
	free(fluxes);
}

void Euler1d::compute_source_term(double** u, double** res, double* Af)
{
	#pragma acc kernels present(u, Af, res)
	{
		#pragma acc loop independent gang worker device_type(nvidia) vector(NVIDIA_VECTOR_LENGTH)
		for(int i = 1; i <= N; i++)
		{
			double p = (g-1.0)*(u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]);
			res[i][1] += p*(Af[i] - Af[i-1]);
		}
	}
}

void Euler1d::apply_boundary_conditions()
{
	double** u = this->u;
	double** prim = this->prim;
	double* bcvalL = this->bcvalL;
	double* bcvalR = this->bcvalR;
	int* bcL = &(this->bcL);
	int* bcR = &(this->bcR);
	#pragma acc parallel present(prim[:N+2][:NVARS], u[:N+2][:NVARS], bcvalL[:NVARS], bcvalR[:NVARS], bcL[:1], bcR[:1]) num_gangs(1)
	{
		if(*bcL == 0)
		{
			// homogeneous Dirichlet - wall
			u[0][0] = u[1][0];
			u[0][1] = -u[1][1];
			u[0][2] = u[1][2];
			prim[0][0] = prim[1][0];
			prim[0][1] = -prim[1][1];
			prim[0][2] = prim[1][2];
		}
		else if(*bcL == 1)
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
				prim[0][0] = rho;
				prim[0][1] = v;
				prim[0][2] = p;
			}
			else if(M_in >= 0)
			{
				// subsonic inflow
				// get conserved variables from pt and Tt specified in bcvalL[0] and bcvalL[1] respectively
				double vold, pold, cold, vold1, pold1, cold1, astar, dpdu, dt0, lambda, du, v, T, p;
				double* uold0 = u[0];
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
				prim[0][0] = u[0][0];
				prim[0][1] = v;
				prim[0][2] = p;

				/*c = sqrt(g*p/u[0][0]);
				M = v/c;
				std::cout << "  apply_boundary_conditions(): Inlet ghost cell mach number = " << M << std::endl;*/
			}
#ifndef _OPENACC
			else
				std::cout << "! Euler1d: apply_boundary_conditions(): Error! Inlet is becoming outlet!" << std::endl;
#endif
		}
		//else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;

		if(*bcR == 0)
		{
			u[N+1][0] = u[N][0];
			u[N+1][1] = -u[N][1];
			u[N+1][2] = u[N][2];
			prim[N+1][0] = prim[N][0];
			prim[N+1][1] = -prim[N][1];
			prim[N+1][2] = prim[N][2];
		}
		else if(*bcR == 3)
		{
			// outflow
			double l1, l2, l3, cold, cold1, pold, pold1, vold, vold1, dt0, r1, r2, r3, Mold, dp, drho, dv, p;
			double* uold = u[N+1];
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
			prim[N+1][0] = u[N+1][0];
			prim[N+1][1] = vold+dv;
			prim[N+1][2] = p;
		}
	//else std::cout << "! Euler1D: apply_boundary_conditions(): BC type not recognized!" << std::endl;
	}
}

Euler1dExplicit::Euler1dExplicit(int num_cells, double length, int leftBCflag, int rightBCflag, std::vector<double> leftBVs, std::vector<double> rightBVs,
		double CFL, std::string inviscid_flux, std::string slope_scheme, std::string face_extrap_scheme, std::string limiter, double f_time, int temporal_order, std::string RKfile)
	: Euler1d(num_cells,length,leftBCflag,rightBCflag,leftBVs,rightBVs, CFL, inviscid_flux, slope_scheme, face_extrap_scheme, limiter), ftime(f_time), temporalOrder(temporal_order)
{
	RKCoeffs = (double**)malloc(temporalOrder*sizeof(double*));
	RKCoeffs[0] = (double*)malloc(temporalOrder*3*sizeof(double));
	for(int i = 0; i < temporalOrder; i++)
		RKCoeffs[i] = *RKCoeffs + i*3;

	std::ifstream rkfile(RKfile);
	for(int i = 0; i < temporalOrder; i++)
		for(int j = 0; j < 3; j++)
			rkfile >> RKCoeffs[i][j];
	rkfile.close();

	std::cout << "Euler1dExplicit: Using " << temporalOrder << "-stage TVD RK scheme; loaded coefficients.\n";
}

Euler1dExplicit::~Euler1dExplicit()
{
	free(RKCoeffs[0]);
	free(RKCoeffs);
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
			prim[i][0] = 1.0;
			prim[i][1] = 0;
			prim[i][2] = 1.0;
		}
		else
		{
			u[i][0] = 0.125;
			u[i][1] = 0;
			u[i][2] = 0.25;
			prim[i][0] = 0.125;
			prim[i][1] = 0;
			prim[i][2] = 0.1;
		}

	double* c = (double*)malloc((N+2)*sizeof(double));
	double** uold = (double**)malloc((N+2)*sizeof(double*));
	//uold[0] = (double*)malloc((N+2)*NVARS*sizeof(double));
	double** ustage = (double**)malloc((N+2)*sizeof(double*));
	ustage[0] = (double*)malloc((N+2)*NVARS*sizeof(double));
	for(int i = 0; i < N+2; i++)
	{
		//uold[i] = *uold + i*NVARS;
		ustage[i] = *ustage + i*NVARS;
	}
	
	for(int i = 0; i < ncell; i++)
		uold[i] = (double*)malloc(NVARS*sizeof(double));

	double** u = this->u;
	double** prim = this->prim;
	double** RKCoeffs = this->RKCoeffs;
	double** dudx = this->dudx;
	double** res = this->res;
	double** prleft = this->prleft;
	double** prright = this->prright;
	double* x = this->x;
	double* dx = this->dx;
	double* A = this->A;
	double* Af = this->Af;
	double* vol = this->vol;
	double* nodes = this->nodes;
	double* bcvalL = this->bcvalL;
	double* bcvalR = this->bcvalR;
	int* bcL = &(this->bcL);
	int* bcR = &(this->bcR);
	//int* N = &(this->N);
	//double* cfl = &(this->cfl);

	#pragma acc enter data copyin(u[:N+2][:NVARS], prim[:N+2][:NVARS], x[:N+2], dx[:N+2], A[:N+2], vol[:N+2], Af[:N+1], nodes[:N+1], RKCoeffs[:temporalOrder][:3], bcvalL[:NVARS], bcvalR[:NVARS], bcL[:1], bcR[:1], g)
	#pragma acc enter data create(dudx[:N+2][:NVARS], res[:N+2][:NVARS], prleft[:N+1][:NVARS], prright[:N+1][:NVARS], dt, mws, istage, c[:N+2], uold[:N+2][:NVARS], ustage[:N+2][:NVARS])

	while(time < ftime)
	{
		std::cout << "Euler1dExplicit: run(): Started time loop" << std::endl;

		#pragma acc update self(u[:N+2][:NVARS])
		std::cout << "Euler1dExplicit: run(): Updated self" << std::endl;

		#pragma acc parallel loop present(u[:N+2][:NVARS], uold[:N+2][:NVARS]) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
		for(int i = 0; i < N+2; i++)
		{
			for(int j = 0; j < NVARS; j++)
			{
				uold[i][j] = u[i][j];
			}
		}
		std::cout << "Euler1dExplicit: run(): Set uold" << std::endl;
		
		// find time step as dt = CFL * min{ dx[i]/(|v[i]|+c[i]) }
		
		mws = 1e10;
		for(int i = 1; i < N+1; i++)
		{
			c[i] = sqrt( g*(g-1.0) * (u[i][2] - 0.5*u[i][1]*u[i][1]/u[i][0]) / u[i][0] );
		}

		//double mws = dx[1]/(fabs(u[1][1]) + c[1]);
		for(int i = 2; i < N+1; i++)
		{
			a = dx[i]/(fabs(u[i][1]) + c[i]);
			if(a < mws) {
				mws = a;
			}
		}

		dt = cfl*mws;

		std::cout << "Euler1dExplicit: run(): Computed dt" << std::endl;

		#pragma acc update device(dt)

		std::cout << "Euler1dExplicit: run(): Updated dt" << std::endl;

		// NOTE: moved apply_boundary_conditions() to the top of the inner loop
		for(istage = 0; istage < temporalOrder; istage++)
		{
			// apply BCs
			{
				apply_boundary_conditions();
			}
			std::cout << "Euler1dExplicit: run():  Applied BCs" << std::endl;

			#pragma acc parallel loop present(u, ustage, res) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
			for(int i = 0; i < N+2; i++)
			{
				for(int j = 0; j < NVARS; j++)
				{
					ustage[i][j] = u[i][j];
					res[i][j] = 0;
				}
			}
			std::cout << "Euler1dExplicit: run():  Set ustage and res" << std::endl;

			cslope->compute_slopes();

			rec->compute_face_values();
			std::cout << "Euler1dExplicit: run():  Computed face values" << std::endl;

			compute_inviscid_fluxes(prleft,prright,res,Af);
			std::cout << "Euler1dExplicit: run():  Computed fluxes" << std::endl;

			compute_source_term(u,res,Af);
			std::cout << "Euler1dExplicit: run():  Computed source terms" << std::endl;

			// RK stage
			#pragma acc parallel loop present(prim, u, uold, ustage, res, vol, dt, RKCoeffs) gang worker vector device_type(nvidia) vector_length(NVIDIA_VECTOR_LENGTH)
			for(int i = 1; i < N+1; i++)
			{
				for(int j = 0; j < NVARS; j++)
					u[i][j] = RKCoeffs[istage][0]*uold[i][j] + RKCoeffs[istage][1]*ustage[i][j] + RKCoeffs[istage][2]*dt/vol[i]*res[i][j];

				prim[i][0] = u[i][0];
				prim[i][1] = u[i][1]/u[i][0];
				prim[i][2] = (g-1.0)*(u[i][2] - 0.5*u[i][1]*prim[i][1]);
			}

		}

		if(step % 10 == 0)
			std::cout << "Euler1dExplicit: run(): Step " << step << " - Time = " << time << std::endl;

		time += dt;
		step++;
	}

	#pragma update self(u[:N+2][:NVARS], prim[:N+2][:NVARS])
	
	#pragma acc exit data delete(u[:N+2][:NVARS], prim[:N+2][:NVARS], x[:N+2], dx[:N+2], A[:N+2], vol[:N+2], Af[:N+1], nodes[:N+1], RKCoeffs[:temporalOrder][:3], bcvalL[:NVARS], bcvalR[:NVARS], bcL, bcR, g, N, cfl, dudx[:N+2][:NVARS], res[:N+2][:NVARS], prleft[:N+1][:NVARS], prright[:N+1][:NVARS], dt, mws, c[:N+2], uold[:N+2][:NVARS], ustage[:N+2][:NVARS])

	free(c);
	for(int i = 0; i < ncell; i++)
		free(uold[i]);
	//free(uold[0]);
	free(uold);
	free(ustage[0]);
	free(ustage);

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
	double cin;

	double pex = bcvalR[0], tex = bcvalR[1], Mex = bcvalR[2], cex, vex;
	
	// set some cells according to inlet condition
	int pn = N/4;
	for(int i = 0; i <= pn; i++)
	{
		u[i][0] = pin/(R*Tin);
		cin = sqrt(g*pin/u[i][0]);
		u[i][1] = u[i][0]*M*cin;
		u[i][2] = pin/(g-1.0)+0.5*u[i][1]*M*cin;
		prim[i][0] = u[i][0];
		prim[i][1] = M*cin;
		prim[i][2] = pin;
	}

	// set last cells according to exit conditions
	for(int i = N+1; i <= N+1; i++)
	{
		u[i][0] = pex/(R*tex);
		cex = sqrt(g*pex/u[i][0]);
		vex = Mex*cex;
		u[i][1] = u[i][0]*vex;
		u[i][2] = pex/(g-1.0) + 0.5*u[i][0]*vex*vex;
		prim[i][0] = u[i][0];
		prim[i][1] = vex;
		prim[i][2] = pex;
	}

	// linearly interpolate cells in the middle
	for(int i = pn; i < N+1; i++)
	{
		for(int j = 0; j < NVARS; j++)
		{
			u[i][j] = (u[N+1][j]-u[pn][j])/(x[N+1]-x[pn])*(x[i]-x[pn]) + u[pn][j];
			prim[i][j] = (prim[N+1][j]-prim[pn][j])/(x[N+1]-x[pn])*(x[i]-x[pn]) + prim[pn][j];
		}
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

		//cslope->compute_slopes();
		rec->compute_face_values();

		compute_inviscid_fluxes(prleft,prright,res,Af);
		compute_source_term(u,res,Af);

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
		{
			for(j = 0; j < NVARS; j++)
				u[i][j] = uold[i][j] + dt[i]/vol[i]*res[i][j];

			prim[i][0] = u[i][0];
			prim[i][1] = u[i][1]/u[i][0];
			prim[i][2] = (g-1.0)*(u[i][2] - 0.5*u[i][1]*prim[i][1]);
		}

		// apply BCs
		apply_boundary_conditions();

		if(step % 10 == 0)
		{
			std::cout << "Euler1dSteadyExplicit: run(): Step " << step << ", relative mass flux norm = " << resnorm/resnorm0 << std::endl;
		}

		step++;
	}

	std::cout << "Euler1dExplicit: run(): Done. Number of time steps = " << step << std::endl;

	/*for(int j = 0; j < NVARS; j++)
	{
		for(int i = 0; i <= N+1; i++)
			std::cout << dudx[i][j] << " ";
		std::cout << std::endl;
	}
	std::cout << std::endl;*/

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
