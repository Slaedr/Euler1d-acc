#include "1deuler.hpp"

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		std::cout << "Euler1d needs a configuration file name. Quitting.\n";
		return -1;
	}

	std::string confile(argv[1]);

	std::vector<double> leftb(NVARS,0);
	std::vector<double> rightb(NVARS,0);
	std::string inv_flux = "llf";
	double cfl = 0.25;
	double f_time = 0.5;
	int temporal_order = 1;

	Euler1dExplicit prob(100, 1.0, 0,0, leftb, rightb, inv_flux, cfl, f_time, temporal_order);

	return 0;
}
