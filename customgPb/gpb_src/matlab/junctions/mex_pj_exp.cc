/*
 * Pj.
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

#include "concurrent/threads/thread.hh"

#include <time.h>
#include <mex.h>


using collections::pointers::auto_collection;
using collections::array_list;
using std::endl;
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

using concurrent::threads::thread;

/********************************** 
 * Matlab matrix conversion routines.
 **********************************/


/*
 * Convert an mxArray to an array.
 */
 
array<double> to_array(const mxArray *a) {
	unsigned long mrows = static_cast<unsigned long>(mxGetM(a));
	unsigned long ncols = static_cast<unsigned long>(mxGetN(a));
	
	if (mrows != 1 && ncols != 1) {
		throw ex_invalid_argument("Array required");
	}
	
	array<double> arr(mrows > ncols ? mrows : ncols);
	double *data = mxGetPr(a);
	for (unsigned long n = 0; n < arr.size(); n++) {
		arr[n] = data[n];
	}
	
	return arr;
}

	typedef array_list<auto_collection< matrix<>, array_list<matrix<> > > > 
		t_2d_matrices;
		
	typedef auto_collection< matrix<>, array_list< matrix<> > >
		t_auto_1d_matrices;
	typedef auto_collection<
		auto_collection< matrix<>, array_list<matrix<> > >,
		array_list<auto_collection< matrix<>, array_list<matrix<> > > > >
		t_auto_2d_matrices;


/* 
 * Convert an mxArray to a matrix.
 */
matrix<> to_matrix(const mxArray *a) {
   unsigned long mrows = static_cast<unsigned long>(mxGetM(a));
   unsigned long ncols = static_cast<unsigned long>(mxGetN(a));
   double *data = mxGetPr(a);
   matrix<> m(mrows, ncols);
   for (unsigned long r = 0; r < mrows; r++) {
      for (unsigned long c = 0; c < ncols; c++) {
         m(r,c) = data[(c*mrows) + r];
      }
   }
   return m;
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

// Shouldn't really be used, convert unsigned long to double
void convertToDouble(const matrix<unsigned long> &in, matrix<double> &out) {
	for (unsigned long r = 0; r < in.size(0); r++) {
		for (unsigned long c = 0; c < in.size(1); c++) {
			out(r, c) = in(r, c);
		}
	}
}


matrix<> &get2dEl(array_list<auto_collection< matrix<>, array_list<matrix<> > > > &coll,
	unsigned long ind1, unsigned long ind2)
{
	return (*(coll[ind1]))[ind2];
}

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


/*
 * Matlab interface.
 * Input arguments are im (in RGB, 3 dim array), gpb, maxo, vect (3D array), params
 * params is an array: labweight, pos_channel_weight, eigvectweight, n_slices, n_oris_gpb, rad_support, min_n_angles, max_n_angles, threshold_rad_support, pos_channel_threshold
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	if (nlhs != 3 || nrhs != 5) {
		cout << "Args!\n";
		return;
	}
	
	//INPUT PARAMETERS
	//should actually come in as rgb, will get converted to Lab real soon
	t_auto_1d_matrices lab_p = to_matrices(prhs[0]);
	t_1d_matrices &lab = *lab_p;
	
	//positive channel and orientations
	matrix<> poschan = to_matrix(prhs[1]);
	matrix<> poschan_ori = to_matrix(prhs[2]);
	
	//parameters for compute_pj
	double *in_params = mxGetPr(prhs[4]);
	if (mxGetNumberOfElements(prhs[4]) != 10) {
		cout << "There should be 10 params" << endl;
		return;
	}
	double labweight = in_params[0];
	double pos_channel_weight = in_params[1];
	double eigvectweight = in_params[2];
	
	unsigned long n_slices = static_cast<unsigned long>(in_params[3]);
	unsigned long n_oris_gpb = static_cast<unsigned long>(in_params[4]);
	unsigned long rad_support = static_cast<unsigned long>(in_params[5]);
	unsigned long min_n_angles = static_cast<unsigned long>(in_params[6]);
	unsigned long max_n_angles = static_cast<unsigned long>(in_params[7]);
	//threshold parameters
	unsigned long threshold_rad_support = static_cast<unsigned long>(in_params[8]);
	double pos_channel_threshold = in_params[9];
	
	//eigvects (do this after params because we need threshold weight)
	t_auto_1d_matrices eigvects;
	unsigned long num_eigvects = 0;
	if (mxGetNumberOfDimensions(prhs[3]) != 3 || eigvectweight <= 1e-9) {
		cout << "ignoring eigvects" << endl;
	} else {
		eigvects = to_matrices(prhs[3]);
		num_eigvects = eigvects->size();
	}
	
	//see if we are going to use lab channel
	unsigned long num_labchan = 3;
	if (labweight <= 1e-9) {
		cout << "ignoring lab" << endl;
		num_labchan = 0;
	}
	
	//compute channel weights from given parameters
	unsigned long total_channels = num_labchan + num_eigvects;
	array<double> channel_weights(total_channels);
	{
		unsigned long curr_channel = 0;
		for (unsigned long i = 0; i < num_labchan; i++) {
			channel_weights[curr_channel] = 1.0 / num_labchan * labweight;
			curr_channel++;
		}
		for (unsigned long i = 0; i < num_eigvects; i++) {
			channel_weights[curr_channel] = 1.0 / num_eigvects * eigvectweight;
			curr_channel++;
		}
	}
	
	//parameters binning & smoothing, not really going to change
	unsigned long num_bins  = 25;                    /* # bins for bg */
	double bg_smooth_sigma = 0.1;
	double cg_smooth_sigma = 0.05;
	double sigma_tg_filt_sm = 2.0;
	double sigma_tg_filt_lg = math::sqrt(2) * 2.0;
	
	unsigned long n_ori = 8;	//for textons TODO different?
	
	matrix<> bg_smooth_kernel = lib_image::gaussian(bg_smooth_sigma * num_bins);
	matrix<> cga_smooth_kernel = lib_image::gaussian(cg_smooth_sigma * num_bins);
	matrix<> cgb_smooth_kernel = lib_image::gaussian(cg_smooth_sigma * num_bins);
	
	//assume everything is same dimension?
	unsigned long im_size_x = poschan.size(0);
	unsigned long im_size_y = poschan.size(1);
	
	//TODO mirror border
	
	// IMAGE PROCESSING
	//convert to grayscale
	matrix<> gray = lib_image::grayscale(lab[0], lab[1], lab[2]);
	//gamma correct
	lib_image::rgb_gamma_correct(lab[0], lab[1], lab[2], 2.5);
	//convert from rgb to Lab
	lib_image::rgb_to_lab(lab[0], lab[1], lab[2]);
	lib_image::lab_normalize(lab[0], lab[1], lab[2]);
	
	//compute filter set
	auto_collection< matrix<>, array_list< matrix<> > > filters_small = 
         lib_image::texton_filters(n_ori, sigma_tg_filt_sm);
    auto_collection< matrix<>, array_list< matrix<> > > filters_large = 
         lib_image::texton_filters(n_ori, sigma_tg_filt_lg);
    array_list< matrix<> > filters;
    filters.add(*filters_small);
    filters.add(*filters_large);
    //compute textons
    auto_collection< matrix<>, array_list< matrix<> > > textons;
    matrix<unsigned long> t_assign = lib_image::textons(
         gray, filters, textons, 64);
	
	// BUILD PARAMETERS FOR PJ
	//accumulate normal channels
	t_auto_1dul_matrices channels_p(new t_1dul_matrices());
	t_1dul_matrices &channels = *channels_p;
	
	auto_ptr<matrix<unsigned long> > toadd_ul;
	
	//L, a, b channels, only add if weight is positive
	for (unsigned long i = 0; i < num_labchan; i++) {
		toadd_ul.reset(new matrix<unsigned long>());
		*toadd_ul = lib_image::quantize_values(lab[i], num_bins);
		channels.add(*toadd_ul);
		toadd_ul.release();
	}
	
	//eigen vector channels, only add if weight is positive
	for (unsigned long i = 0; i < num_eigvects; i++) {
		toadd_ul.reset(new matrix<unsigned long>());
		*toadd_ul = lib_image::quantize_values((*eigvects)[i], num_bins);
		channels.add(*toadd_ul);
		toadd_ul.release();
	}
	
	//outputs... TODO assume all are same
	matrix<> pjs(im_size_x, im_size_y);
	matrix< array<double> > jangles(im_size_x, im_size_y);
	matrix<> init_ests(im_size_x, im_size_y);
	
	try {
		lib_image::compute_pj_exp(
			channels, poschan, poschan_ori,
			channel_weights, pos_channel_weight,
			threshold_rad_support, pos_channel_threshold,
			n_slices, n_oris_gpb, rad_support, min_n_angles, max_n_angles, num_bins,
			pjs, jangles, init_ests,
			vector(bg_smooth_kernel));
	} catch (exception e) {
		cout << e;
		e.raise();
	}
	
	//CONVERT JANGLES to CELL ARRAY
	//jangles is a matrix of different sized arrays
	mxArray *cangles = mxCreateCellMatrix(
		static_cast<int>(im_size_x),
		static_cast<int>(im_size_y));
	unsigned long i = 0;
	for (unsigned long y = 0; y < im_size_y; y++) {
	for (unsigned long x = 0; x < im_size_x; x++) {
			int n_angles = static_cast<int>(jangles(x, y).size());
			mxArray *a = mxCreateDoubleMatrix(1, n_angles, mxREAL);
			double *data = mxGetPr(a);
			for (int j = 0; j < n_angles; j++) {
				data[j] = jangles(x, y)[j];
			}
			mxSetCell(cangles, static_cast<int>(i), a);
			i++;
		}
	} 
	
	//return stuff
	plhs[0] = to_mxArray(pjs);
	plhs[1] = cangles;
	plhs[2] = to_mxArray(init_ests);
}










