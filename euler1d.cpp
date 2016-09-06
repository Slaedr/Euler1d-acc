#include "1deuler.hpp"

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		std::cout << "Euler1d needs a configuration file name. Quitting.\n";
		return -1;
	}

	std::string confile(argv[1]);
	std::ifstream conf(confile);

	int leftbc, rightbc, temporal_order, N, areatype;
	std::vector<double> leftbv(NVARS,0);
	std::vector<double> rightbv(NVARS,0);
	std::string inv_flux, areafile, outputfile, dum;
	double cfl, f_time, L;

	conf >> dum; conf >> outputfile;
	conf >> dum; conf >> N;
	conf >> dum; conf >> L;
	conf >> dum; conf >> leftbc;
	conf >> dum; conf >> rightbc;
	conf >> dum;
	for(int i = 0; i < NVARS; i++)
		conf >> leftbv[i];
	conf >> dum;
	for(int i = 0; i < NVARS; i++)
		conf >> rightbv[i];
	conf >> dum; conf >> inv_flux;
	conf >> dum; conf >> cfl;
	conf >> dum; conf >> f_time;
	conf >> dum; conf >> temporal_order;
	conf >> dum; conf >> areatype;
	if(areatype != 0)
	{
		conf >> dum; conf >> areafile;
	}
	conf.close();
	
	std::vector<double> plist, areas(1,1.0);

	Euler1dExplicit prob(N, L, leftbc, rightbc, leftbv, rightbv, inv_flux, cfl, f_time, temporal_order);
	prob.generate_mesh(0,plist);
	prob.set_area(0,areas);
	std::cout << N << " " << L << " " << leftbc << " " << rightbc << " " << inv_flux << " " << cfl << " " << f_time << " " << temporal_order << std::endl;

	prob.run();

	prob.postprocess(outputfile);

	return 0;
}
