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
class SlopeLimiter1
{
protected:
	const double eps;							///< Small number needed by some limiters
public:
	SlopeLimiter1(double _eps);
	virtual double limiter_function(double a, double b) const = 0;
};

/// `Limiter' with value 1.0 always
class NoLimiter1 : public SlopeLimiter1
{
public:
	NoLimiter1(double eps);
	double limiter_function(double a, double b) const;
};

/// Minmod limiter
class MinmodLimiter1 : public SlopeLimiter1
{
public:
	MinmodLimiter1(double eps);
	double limiter_function(double a, double b) const;
};

/// Van Albada limiter
class VanAlbadaLimiter1 : public SlopeLimiter1
{
public:
	VanAlbadaLimiter1(double eps);
	double limiter_function(double a, double b) const;
};


/// Interface for limiters which take a single argument - ratio of consecutive differences
class Limiter
{
public:
	virtual double limiter_function(double r) const = 0;
};

/// For no limiting
class NoLimiter : public Limiter
{
public:
	double limiter_function(double r) const;
};

/// Val Albada (alternative form) limiter
/** Phi(r) := 2r / (1+r^2).
 */
class VanAlbadaLimiter : public Limiter
{
public:
	double limiter_function(double r) const;
};

inline double VanAlbadaLimiter::limiter_function(double r) const
{
	return 2.0*r/(1.0+r*r);
}

/// Limiter by Hemker and Koren; see Blazek's book.
/** Phi(r) := 3r / (2r^2 - r + 2).
 */
class HemkerKorenLimiter : public Limiter
{
public:
	double limiter_function(double r) const;
};

class MinmodLimiter : public Limiter
{
public:
	double limiter_function(double r) const;
};

/// Van Leer limiter - bad
class VanLeerLimiter : public Limiter
{
public:
	double limiter_function(double r) const;
};

#pragma acc routine seq
double vanalbada_limiter_function(double r);

#endif
