/** \file limiters.cpp
 * \brief Implementation of slope limiters
 */

#include "limiters.hpp"

SlopeLimiter1::SlopeLimiter1(double _eps) : eps(_eps)
{ }

NoLimiter1::NoLimiter1(double eps) : SlopeLimiter1(eps)
{ }

double NoLimiter1::limiter_function(double a, double b) const
{
	return a<b ? a : b;
}

MinmodLimiter1::MinmodLimiter1(double eps) : SlopeLimiter1(eps)
{ }

double MinmodLimiter1::limiter_function(double a, double b) const
{
	if(fabs(a) <= fabs(b) && a*b >= 0)
		return a;
	else if(fabs(a) > fabs(b) && a*b >= 0)
		return b;
	else
		return 0.0;
}

VanAlbadaLimiter1::VanAlbadaLimiter1(double eps) : SlopeLimiter1(eps)
{ }

double VanAlbadaLimiter1::limiter_function(double a, double b) const
{
	double val = (2.0*a*b + eps)/(a*a + b*b + eps);
	if(val > 0.0) return val;
	else return 0.0;
}


double NoLimiter::limiter_function(double r) const
{
	return 1.0;
}

double VanAlbadaLimiter::limiter_function(double r) const
{
	return 2.0*r/(1.0+r*r);
}

double HemkerKorenLimiter::limiter_function(double r) const
{
	return 3.0*r/(2.0*r*r-r+2.0);
}
