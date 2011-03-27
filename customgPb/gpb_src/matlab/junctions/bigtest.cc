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



matrix<> weight_matrix_disc(unsigned long r, long r_inner = -1) {
   /* initialize matrix */
   unsigned long size = 2*r + 1;
   matrix<> weights(size, size);
   /* set values in disc to 1 */
   long radius = static_cast<long>(r);
   long r_sq = radius * radius;
   long r_inner_sq = r_inner < 0 ? r_inner : r_inner * r_inner;
   unsigned long ind = 0;
   for (long x = -radius; x <= radius; x++) {
      long x_sq = x * x;
      for (long y = -radius; y <= radius; y++) {
         /* check if index is within disc */
         long y_sq = y * y;
         if ((x_sq + y_sq) <= r_sq && (x_sq + y_sq) > r_inner_sq)
            weights[ind] = 1;
         /* increment linear index */
         ind++;
      }
   }
   
   return weights;
}

void testWeights() {
	cout << weight_matrix_disc(5) << endl;
	cout << weight_matrix_disc(5, 0) << endl;
	cout << weight_matrix_disc(5, 0.5) << endl;
	cout << weight_matrix_disc(5, 1) << endl;
	cout << weight_matrix_disc(5, 2) << endl;
}


/*
 * Construct orientation slice lookup map, have slices go all the way
 * around instead of mirroring across origin
 */
matrix<unsigned long> half_orientation_slice_map(
   unsigned long size_x, unsigned long size_y, unsigned long n_ori)
{
   /* initialize map */
   matrix<unsigned long> slice_map(size_x, size_y);
   /* compute orientation of each element from center */
   unsigned long ind = 0;
   double x = -static_cast<double>(size_x - 1)/2;
   for (unsigned long n_x = 0; n_x < size_x; n_x++) {
      double y = -static_cast<double>(size_y - 1)/2;
      for (unsigned long n_y = 0; n_y < size_y; n_y++) {
         /* compute orientation index */
         double ori = std::atan2(y, x) + M_PIl;
         long idx = static_cast<long>(
            math::ceil(ori / (2 * M_PIl) * static_cast<double>(n_ori))
         ) - 1;
         
         slice_map[ind] = idx;
         /* increment position */
         ind++;
         y++;
      }
      /* increment x-coordinate */
      x++;
   }
   return slice_map;
}

matrix<unsigned long> half_orientation_slice_map_offset(
	unsigned long size_x, unsigned long size_y, unsigned long n_ori)
{
	matrix<unsigned long> slice_map =
		half_orientation_slice_map(size_x, size_y, n_ori * 2);
	for (unsigned long x = 0; x < size_x; x++) {
	for (unsigned long y = 0; y < size_y; y++) {
		if (slice_map(x, y) % 2 == 0)
			slice_map(x, y) /= 2;
		else if (slice_map(x, y) == 2 * n_ori - 1)
			slice_map(x, y) = 0;
		else
			slice_map(x, y) = (slice_map(x, y) + 1) / 2;
	}
	}
	
	return slice_map;
}

void testSlices() {
	unsigned long n_slices = 8;
	unsigned long rad = 3;
	
	array<unsigned long> counts(n_slices), counts2(n_slices);
	matrix<> w = weight_matrix_disc(rad, true);
	matrix<unsigned long> sm = half_orientation_slice_map(w.size(0), w.size(1), n_slices);
	matrix<unsigned long> sm2 = half_orientation_slice_map_offset(w.size(0), w.size(1), n_slices);
	for (unsigned long x = 0; x < sm.size(0); x++) {
	for (unsigned long y = 0; y < sm.size(1); y++) {
		if (w(x, y) == 0) sm(x, y) = 0;
		else counts[sm(x, y)]++;
		
		if (w(x, y) == 0) sm2(x, y) = 0;
		else counts2[sm2(x, y)]++;
	}
	}
	cout << sm << endl
		<< counts << endl
		<< sm2 << endl
		<< counts2 << endl;
}

unsigned long linear_index(
   const array<unsigned long>& dims, 
   unsigned long x, 
   unsigned long y,
   unsigned long z)
{
   unsigned long n_dims = dims.size();
   unsigned long y_size = (n_dims > 1) ? dims[1] : 1;
   unsigned long z_size = (n_dims > 2) ? dims[2] : 1;
   unsigned long index = (x*y_size + y) * z_size + z;
   for (unsigned long n = 3; n < n_dims; n++)
      index *= dims[n];
   return index;
}

unsigned long linear_index(
   const array<unsigned long>& dims,
   const array<unsigned long>& i)
{
   unsigned long n_dims     = dims.size();
   unsigned long n_dims_i   = i.size();
   unsigned long n_dims_min = (n_dims < n_dims_i) ? n_dims : n_dims_i;
   unsigned long index = 0;
   for (unsigned long n = 0; n < n_dims_min; n++)
      index = index*dims[n] + i[n];
   for (unsigned long n = n_dims_min; n < n_dims; n++)
      index *= dims[n];
   for (unsigned long n = n_dims_min; n < n_dims_i; n++)
      index += i[n];
   return index;
}

typedef auto_collection< matrix<>, array_list< matrix<> > >
	t_auto_1d_matrices;
typedef auto_collection<
	auto_collection< matrix<>, array_list<matrix<> > >,
	array_list<auto_collection< matrix<>, array_list<matrix<> > > > >
	t_auto_2d_matrices;
	
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
}

void templateFor2dMatrices() {

/*
	auto_collection<
		auto_collection< matrix<>, array_list<matrix<> > >,
		array_list<auto_collection< matrix<>, array_list<matrix<> > > > > coll(
		new array_list<auto_collection< matrix<>, array_list<matrix<> > > >());
		
	for (int i = 0; i < 5; i++) {
		auto_ptr<auto_collection< matrix<>, array_list< matrix<> > > > ptr(
			new auto_collection< matrix<>, array_list< matrix<> > > (
				new array_list< matrix<> > ()));
		for (int j = 0; j < 8; j++) {
			auto_ptr<matrix<> > iptr(new matrix<>());
			(**ptr).add(*iptr);
			iptr.release();
		}
		coll->add(*ptr);
		ptr.release();
	}
	
	matrix<> &m = (*(*coll)[2])[3]; */
}

typedef array_list<auto_collection< matrix<>, array_list<matrix<> > > > 
	t_2d_matrices;
	
typedef auto_collection< matrix<>, array_list< matrix<> > >
	t_auto_1d_matrices;
typedef auto_collection<
	auto_collection< matrix<>, array_list<matrix<> > >,
	array_list<auto_collection< matrix<>, array_list<matrix<> > > > >
	t_auto_2d_matrices;


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

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	testWeights();
}
