/** \file limiters.cpp
 * \brief Implementation of slope limiters
 */

#include "limiters.hpp"

SlopeLimiter::SlopeLimiter(double _eps) : eps(_eps)
{ }

NoLimiter::NoLimiter(double eps) : SlopeLimiter(eps)
{ }

double NoLimiter::limiter_function(double a, double b) const
{
	return 1.0;
}

VanAlbadaLimiter::VanAlbadaLimiter(double eps) : SlopeLimiter(eps)
{ }

double VanAlbadaLimiter::limiter_function(double a, double b) const
{
	double val = (2.0*a*b + eps)/(a*a + b*b + eps);
	if(val > 0.0) return val;
	else return 0.0;
}
