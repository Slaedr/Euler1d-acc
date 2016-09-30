/** \file limiters.hpp
 * \brief Limiters for TVD/MUSCL schemes
 * \author Aditya Kashi
 * \date 20 September 2016
 */

#ifndef __LIMITERS_H

#define __LIMITERS_H

#ifndef __DEFINITIONS_H
#include "definitions.hpp"
#endif

/// Interface for TVD/MUSCL slope limiters
class SlopeLimiter
{
protected:
	const double eps;							///< Small number needed by some limiters
public:
	SlopeLimiter(double _eps);
	virtual double limiter_function(double a, double b) const = 0;
};

/// `Limiter' with value 1.0 always
class NoLimiter : public SlopeLimiter
{
public:
	NoLimiter(double eps);
	double limiter_function(double a, double b) const;
};

/// Minmod limiter
class MinmodLimiter : public SlopeLimiter
{
public:
	MinmodLimiter(double eps);
	double limiter_function(double a, double b) const;
};

/// Van Albada limiter
class VanAlbadaLimiter : public SlopeLimiter
{
public:
	VanAlbadaLimiter(double eps);
	double limiter_function(double a, double b) const;
};

#endif
