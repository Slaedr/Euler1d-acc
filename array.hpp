#ifndef __ARRAY_H
#define __ARRAY_H

#include <cstdio>

template<typename T>
class Array1d
{
protected:
	T* data;
	size_t size;
	bool alloc;
public:
	Array1d(size_t sz)
	{
		size = sz;
		data = new T[size];
		alloc = true;
	}
	~Array1d()
	{
		if(alloc)
			delete [] data;
	}

	T& operator[](const size_t i)
	{
		return data[i];
	}

	T operator[](const size_t i) const
	{
		return data[i];
	}
};

/** \brief Row-major 2D array storage
 *
 * Stores members in a contiguous 1D block.
 */
template <typename T>
class Array2d
{
protected:
	T* data;
	size_t nrows;
	size_t ncols;
	size_t size;
	bool alloc;
public:
	Array2d(size_t nr, size_t nc)
	{
		nrows = nr; ncols = nc;
		size = ncols*nrows;
		data = new T[size];
		alloc = true;
	}
	~Array2d()
	{
		if(alloc)
			delete [] data;
	}

	T& operator[][](const size_t i, const size_t j)
	{
		return data[i*nrows+j];
	}

	T operator[][](const size_t i, const size_t j) const
	{
		return data[i*nrows+j];
	}
};

#endif
