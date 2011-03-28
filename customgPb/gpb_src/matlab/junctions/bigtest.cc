/*
 * Pb. 
 */
#include "collections/pointers/auto_collection.hh"
#include "collections/array_list.hh"
#include "io/streams/cout.hh"
#include "io/streams/iomanip.hh"
#include "io/streams/ios.hh"
#include "lang/array.hh"
#include "lang/exceptions/exception.hh"
#include "lang/exceptions/ex_index_out_of_bounds.hh"
#include "lang/exceptions/ex_invalid_argument.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/libraries/lib_image.hh"
#include "math/math.hh"
#include "math/matrices/matrix.hh"
#include "math/matrices/functors/matrix_dispersion_functors.hh"

#include <time.h>
#include <mex.h>


using collections::pointers::auto_collection;
using collections::array_list;
using io::streams::cout;
using io::streams::ios;
using io::streams::iomanip::setiosflags;
using io::streams::iomanip::setw;
using lang::array;
using lang::exceptions::exception;
using lang::exceptions::ex_index_out_of_bounds;
using lang::exceptions::ex_invalid_argument;
using lang::pointers::auto_ptr;
using math::libraries::lib_image;
using math::matrices::matrix;
using math::matrices::vector;
using math::matrices::functors::dispersion_functors;

using std::cout;
using std::endl;


/*

void testDisp() {

	t_auto_1d_matrices slice_hists_p(new array_list<matrix<> >());
	array_list<matrix<> > slice_hists;
	for (unsigned long i_angle = 0; i_angle < 10; i_angle++) {
		auto_ptr<matrix<> > angle_ptr(new matrix<>(1, 3));
		slice_hists.add(*angle_ptr);
		angle_ptr.release();
	}
	slice_hists[0](0, 0) = 1; 
	slice_hists[1](0, 0) = 1;
	slice_hists[2](0, 0) = 1;
	slice_hists[3](0, 0) = 1;
	slice_hists[4](0, 0) = 0;
	slice_hists[5](0, 0) = 0;
	slice_hists[6](0, 0) = 0;
	slice_hists[7](0, 0) = 0;
	slice_hists[8](0, 0) = 0;
	slice_hists[9](0, 0) = 0;
	
	slice_hists[0](0, 1) = 0; 
	slice_hists[1](0, 1) = 0;
	slice_hists[2](0, 1) = 0;
	slice_hists[3](0, 1) = 0;
	slice_hists[4](0, 1) = 1;
	slice_hists[5](0, 1) = 1;
	slice_hists[6](0, 1) = 1;
	slice_hists[7](0, 1) = 1;
	slice_hists[8](0, 1) = 1;
	slice_hists[9](0, 1) = 1; 
	
	slice_hists[0](0, 2) = 0; 
	slice_hists[1](0, 2) = 0;
	slice_hists[2](0, 2) = 0;
	slice_hists[3](0, 2) = 0;
	slice_hists[4](0, 2) = 1;
	slice_hists[5](0, 2) = 1;
	slice_hists[6](0, 2) = 1;
	slice_hists[7](0, 2) = 1;
	slice_hists[8](0, 2) = 1;
	slice_hists[9](0, 2) = 1; 
	
	
	array_list<matrix<> > sub;
	slice_hists.subarray(4, 7, sub);

	for (unsigned long i = 0; i < sub.size(); i++) {
		cout << sub[i] << endl;
	}
	cout << dispersion_functors<double>::coll_matrix_var_L2()(sub) << endl;
} */


typedef array_list<matrix<> > t_1d_matrices;
typedef auto_collection<matrix<>, t_1d_matrices> t_auto_1d_matrices;
typedef array_list<t_auto_1d_matrices> t_2d_matrices;
typedef auto_collection<t_auto_1d_matrices,	t_2d_matrices> t_auto_2d_matrices;

typedef array_list<matrix<unsigned long> > t_1dul_matrices;
typedef auto_collection<matrix<unsigned long>, t_1dul_matrices> t_auto_1dul_matrices;
	

t_auto_1d_matrices to_matrices(const mxArray *a) {
	//assume that
	if (mxGetNumberOfDimensions(a) != 3) {
		throw ex_invalid_argument("Number of dims not 3");
	}
	
	const mwSize *dims = mxGetDimensions(a);
	unsigned long mrows = dims[0];
	unsigned long ncols = dims[1];
	unsigned long pmats = dims[2];
	double *data = mxGetPr(a);

	t_auto_1d_matrices out(new array_list<matrix<> >());
	for (unsigned long k = 0; k < pmats; k++) {
		auto_ptr<matrix<> > mptr(new matrix<>(mrows, ncols));
		
		
		for (unsigned long r = 0; r < mrows; r++) {
			for (unsigned long c = 0; c < ncols; c++) {
				(*mptr)(r, c) = data[(k*ncols + c)*mrows + r];
			}
		}
		
		out->add(*mptr);
		mptr.release();
	}
	
	return out;
}


// Shouldn't really be used, convert unsigned long to double
void convertToDouble(const matrix<unsigned long> &in, matrix<double> &out) {
	for (unsigned long r = 0; r < in.size(0); r++) {
		for (unsigned long c = 0; c < in.size(1); c++) {
			out(r, c) = in(r, c);
		}
	}
}



/*
 * Convert a 2D matrix to an mxArray.
 */
mxArray* to_mxArray(const matrix<>& m) {
   unsigned long mrows = m.size(0);
   unsigned long ncols = m.size(1);
   mxArray *a = mxCreateDoubleMatrix(
      static_cast<int>(mrows),
      static_cast<int>(ncols),
      mxREAL
   );
   double *data = mxGetPr(a);
   for (unsigned long r = 0; r < mrows; r++) {
      for (unsigned long c = 0; c < ncols; c++) {
         data[(c*mrows) + r] = m(r,c);
      }
   }
   return a;
}

	double var(const matrix<double> &m) {
		/* matlab code, does identical computation (well, stdev)
		m = sum(ushist .* (0:(nbins - 1)));
		sd = sqrt(sum(ushist' .* ((0:(nbins - 1))' - m) .^ 2));
		m = m / (nbins - 1);
		sd = sd / (nbins - 1);
		*/
		
		//compute mean first
		double mean = 0;
		double variance = 0;
		unsigned long dim = m.size(0);
		
		for (unsigned long i = 0; i < dim; i++) {
			mean += static_cast<double>(m(i, 0)) * (static_cast<double>(i) + 0.5);
		}
		
		//now variance
		for (unsigned long i = 0; i < dim; i++) {
			double diff = (static_cast<double>(i) + 0.5) - mean;
			variance += static_cast<double>(m(i, 0)) * diff * diff;
		}
		cout << mean << ", " << variance << endl;
		
		variance /= (dim - 1) * (dim - 1);
		
		return variance;
	}

void testVar() {
	matrix<double> x(25, 1);
	x(0, 0) = 0.1636;
	x(1, 0) = 0.1394;
	x(2, 0) = 0.1013;
	x(3, 0) = 0.0627;
	x(4, 0) = 0.0330;
	x(5, 0) = 0;
	x(6, 0) = 0;
	x(7, 0) = 0;
	x(8, 0) = 0;
	x(9, 0) = 0;
	x(10, 0) = 0;
	x(11, 0) = 0;
	x(12, 0) = 0;
	x(13, 0) = 0;
	x(14, 0) = 0;
	x(15, 0) = 0;
	x(16, 0) = 0;
	x(17, 0) = 0;
	x(18, 0) = 0;
	x(19, 0) = 0;
	x(20, 0) = 0.0330;
	x(21, 0) = 0.0627;
	x(22, 0) = 0.1013;
	x(23, 0) = 0.1394;
	x(24, 0) = 0.1636;
	
	double sum = 0;
	for (unsigned long i = 0; i < 25; i++) {
		sum += x(i, 0);
	}
	for (unsigned long i = 0; i < 25; i++) {
		x(i, 0) /= sum;
	}
	
	
	matrix<double> y(5, 1);
	y(0, 0) = 1;
	y(1, 0) = 0;
	y(2, 0) = 0;
	y(3, 0) = 0;
	y(4, 0) = 1;
	
	sum = 0;
	for (unsigned long i = 0; i < 5; i++) {
		sum += y(i, 0);
	}
	for (unsigned long i = 0; i < 5; i++) {
		y(i, 0) /= sum;
	}
	
	cout << x << endl << y << endl;
	
	cout << var(x) << endl;
	cout << var(y) << endl;
	
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	testVar();
	
}
