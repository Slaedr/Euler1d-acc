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

	int leftbc, rightbc, temporal_order, N;
	std::vector<double> leftbv(NVARS,0);
	std::vector<double> rightbv(NVARS,0);
	std::string inv_flux, dum;
	double cfl, f_time, L;

	conf >> dum; conf >> N;
	conf >> dum; conf >> L;
	conf >> dum; conf >> leftbc;
	conf >> dum; conf >> rightbc;
	conf dum;
	for(int i = 0; i < NVARS; i++)
		conf >> leftbv[i];
	conf >> dum;
	for(int i = 0; i < NVARS; i++)
		conf >> rightbv[i];
	conf >> dum; conf >> inv_flux;
	conf >> dum; conf >> cfl;
	conf >> dum; conf >> f_time;
	conf >> dum; conf >> temporal_order;
	conf.close();

	Euler1dExplicit prob(N, L, leftbc, rightbc, leftbv, rightbv, inv_flux, cfl, f_time, temporal_order);

	return 0;
}
