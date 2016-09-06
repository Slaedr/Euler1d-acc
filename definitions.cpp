#include "definitions.hpp"

void matprint(const std::vector<std::vector<double>>& mat)
{
	std::cout << std::endl;
	for(int i = 0; i < mat.size(); i++)
	{
		for(int j = 0; j < mat[i].size(); j++)
			std::cout << mat[i][j] << " \t";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

