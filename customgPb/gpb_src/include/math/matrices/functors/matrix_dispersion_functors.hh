/*
 * Functors for dispersion statistics
 */
 
#ifndef MATH__MATRICES__FUNCTORS__MATRIX_DISPERSION_FUNCTORS_HH
#define MATH__MATRICES__FUNCTORS__MATRIX_DISPERSION_FUNCTORS_HH

#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "math/matrices/functors/matrix_distance_functors.hh"
#include "functors/distanceable_functors.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/iterators/iterator.hh"
#include "lang/pointers/auto_ptr.hh"
#include "collections/abstract/collection.hh"

namespace math {
namespace matrices {
namespace functors {

using functors::distanceable_functor;
using lang::exceptions::ex_invalid_argument;
using lang::iterators::iterator;
using lang::pointers::auto_ptr;
using math::matrices::functors::matrix_distance_functors;
using math::matrices::matrix;
using collections::abstract::collection;



////dispersion functors on collections of histograms
template <typename T>
class dispersion_functor {
public:
	virtual ~dispersion_functor() { }
	virtual double operator()(const T &) const = 0;
};


template <typename T>
class dispersion_functor<const T> : public dispersion_functor<T> { };


template <typename T>
class collection_matrix_var : public dispersion_functor<const collection<matrix<T> > > {
protected:
	const distanceable_functor<matrix<double>, double> &_f_norm;
public:
	collection_matrix_var(const distanceable_functor<matrix<double>, double> &f_norm)
		: _f_norm(f_norm)
		{}

	double operator()(const collection<matrix<T> > &ms) const {
		//checks
		auto_ptr<iterator<matrix<T> > > i = ms.iter_create();
		
		if (!i->has_next()) {
			throw ex_invalid_argument("Empty collection in collection_matrix_var");
		}
		//else
		matrix<T> &first = i->next();
//		array<unsigned long> dims = first.dimensions();
		unsigned long dim1 = first.size(0);
		unsigned long dim2 = first.size(1);
		//assume dims is non-zero and all dims are same
		
		//compute mean
		matrix<double> mean(dim1, dim2);
		matrix<double> zeros(dim1, dim2);
		//fill this in:
		for (unsigned long n = 0; n < first._size; n++)
			mean._data[n] = static_cast<double>(first._data[n]);
		unsigned long count = 1;
		
		while (i->has_next()) {
			count++;
			
			//do addition by hand because we need to typecast
			matrix<T> &curr = i->next();
			for (unsigned long n = 0; n < curr._size; n++) {
				mean._data[n] += static_cast<double>(curr._data[n]);
			}
		}
		mean /= static_cast<double>(count);
		
		//variance
		double var = 0;
		matrix<double> diff(dim1, dim2);
		
		i.reset();
		i = ms.iter_create();
		
		while (i->has_next()) {
			matrix<T> &curr = i->next();
			for (unsigned long n = 0; n < curr._size; n++) {
				diff._data[n] = mean._data[n] - static_cast<double>(curr._data[n]);
			}
			var += _f_norm(diff, zeros);
		}
		
		var /= static_cast<double>(count);
		return var;
	}
};



/////dispersion functors on a single histogram
template <typename T>
class hist_disp_functor {
public:
	virtual ~hist_disp_functor() { }
	virtual double operator()(const T &) const = 0;
};

template <typename T>
class hist_disp_functor<const T> : public hist_disp_functor<T> { };


//interprets the matrix as a one-dimensional histogram going from 0 to 1
template <typename T>
class hist_matrix_var : public hist_disp_functor<const matrix<T> > {
public:
	hist_matrix_var() {}
	
	double operator()(const matrix<T> &m) const {
		/* matlab code, does identical computation (well, stdev)
		m = sum(ushist .* (0:(nbins - 1)));
		sd = sqrt(sum(ushist' .* ((0:(nbins - 1))' - m) .^ 2));
		m = m / (nbins - 1);
		sd = sd / (nbins - 1);
		*/
		
		//compute mean first
		double mean = 0;
		double variance = 0;
		unsigned long dim = m._size;
		
		for (unsigned long i = 0; i < dim; i++) {
			mean += static_cast<double>(m._data[i]) * static_cast<double>(i);
		}
		
		//now variance
		for (unsigned long i = 0; i < dim; i++) {
			double diff = static_cast<double>(i) - mean;
			variance += static_cast<double>(m._data[i]) * diff * diff;
		}
		variance /= (dim - 1) * (dim - 1);
		
		return variance;
	}
};




//globally accessible dispersion functors
template <typename T>
class dispersion_functors {
public:
	static const collection_matrix_var<T> &coll_matrix_var_L2();
	static const hist_matrix_var<T> &hist_var();
};

template <typename T>
const collection_matrix_var<T> &dispersion_functors<T>::coll_matrix_var_L2() {
	static const collection_matrix_var<T> *f =
		new collection_matrix_var<T>(matrix_distance_functors<double>::L2_distance());
	return *f;
}

template <typename T>
const hist_matrix_var<T> &dispersion_functors<T>::hist_var() {
	static const hist_matrix_var<T> *f = new hist_matrix_var<T>();
	return *f;
}



} //functors
} //matrices
} //math

#endif
