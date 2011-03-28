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
	
	
	cout << "hey there" << endl;
	double x;
	std::cin >> x;
	
	//INPUT PARAMETERS
	//should actually come in as rgb, will get converted to Lab real soon
	t_auto_1d_matrices lab_p = to_matrices(prhs[0]);
	t_1d_matrices &lab = *lab_p;
	
	//positive channel and orientations
	matrix<> poschan = to_matrix(prhs[1]);
	matrix<> poschan_ori = to_matrix(prhs[2]);
	
	//parameters for compute_pj
	double *in_params = mxGetPr(prhs[4]);
	if (mxGetNumberOfElements(prhs[4]) != 18) {
		cout << "There should be 18 params" << endl;
		return;
	}
	double l_weight = in_params[0];
	double ab_weight = in_params[1];
	double pos_channel_weight = in_params[2];
	double eigvect_weight = in_params[3];
	double text_weight = in_params[4];
	
	double l_homog_weight = in_params[5];
	double ab_homog_weight = in_params[6];
	double eigvect_homog_weight = in_params[7];
	double text_homog_weight = in_params[8];
	
	unsigned long n_slices = static_cast<unsigned long>(in_params[9]);
	unsigned long n_oris_gpb = static_cast<unsigned long>(in_params[10]);
	unsigned long rad_support = static_cast<unsigned long>(in_params[11]);
	double rad_inner_support = in_params[12];
	unsigned long min_n_angles = static_cast<unsigned long>(in_params[13]);
	unsigned long max_n_angles = static_cast<unsigned long>(in_params[14]);
	//threshold parameters
	unsigned long threshold_rad_support = static_cast<unsigned long>(in_params[15]);
	double pos_channel_threshold = in_params[16];
	double norm_power = in_params[17];
	
	
	
	
	//STATIC PARAMETERS
	//parameters for binning & smoothing, not really going to change
	unsigned long num_bins  = 25;                    /* # bins for bg */
	double bg_smooth_sigma = 0.1;
	double cg_smooth_sigma = 0.05;
	double sigma_tg_filt_sm = 2.0;
	double sigma_tg_filt_lg = math::sqrt(2) * 2.0;
	unsigned long num_textons = num_bins; //used to be 64.... but would require big changes
	
	unsigned long n_ori = 8;	//for textons TODO different?
	
	matrix<> bg_smooth_kernel = lib_image::gaussian(bg_smooth_sigma * num_bins);
	matrix<> cga_smooth_kernel = lib_image::gaussian(cg_smooth_sigma * num_bins);
	matrix<> cgb_smooth_kernel = lib_image::gaussian(cg_smooth_sigma * num_bins);
	
	//assume everything is same dimension?
	unsigned long im_size_x = poschan.size(0);
	unsigned long im_size_y = poschan.size(1);
	
	
	
	// IMAGE PROCESSING
	//convert to grayscale
	matrix<> gray = lib_image::grayscale(lab[0], lab[1], lab[2]);
	//gamma correct
	lib_image::rgb_gamma_correct(lab[0], lab[1], lab[2], 2.5);
	//convert from rgb to Lab
	lib_image::rgb_to_lab(lab[0], lab[1], lab[2]);
	lib_image::lab_normalize(lab[0], lab[1], lab[2]);
	
	//for texture channel
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
         gray, filters, textons, num_textons);
         
         
	
    //GET EIGENVECTORS
	//eigvects (do this after params because we need threshold weight)
	t_auto_1d_matrices eigvects;
	unsigned long num_eigvects = 0;
	if (mxGetNumberOfDimensions(prhs[3]) != 3) {
		cout << "eigvects incorrect" << endl;
		eigvect_weight = eigvect_homog_weight = 0;
	} else {
		eigvects = to_matrices(prhs[3]);
		num_eigvects = eigvects->size();
	}
         
	//COMPILE CHANNELS AND WEIGHTS
	//channels
	t_auto_1dul_matrices channels_p(new t_1dul_matrices());
	t_1dul_matrices &channels = *channels_p;
	
	//channel weights
	array<double> channel_weights(num_eigvects + 4);
	array<double> homog_weights(num_eigvects + 4);

	unsigned long total_channels = 0;	
	
	auto_ptr<matrix<unsigned long> > toadd_ul; //temp variable for converting
	
	//L, a, b channels, only add if weight is positive
	if (l_weight < 1e-9 && l_homog_weight < 1e-9) {
		cout << "Ignoring L" << endl;
	} else {
		toadd_ul.reset(new matrix<unsigned long>());
		*toadd_ul = lib_image::quantize_values(lab[0], num_bins);
		channels.add(*toadd_ul);
		toadd_ul.release();
		
		channel_weights(total_channels) = l_weight;
		homog_weights(total_channels) = l_homog_weight;
		
		total_channels++;
	}
	
	if (ab_weight < 1e-9 && ab_homog_weight < 1e-9) {
		cout << "Ignoring ab" << endl;
	} else {
		toadd_ul.reset(new matrix<unsigned long>());
		*toadd_ul = lib_image::quantize_values(lab[1], num_bins);
		channels.add(*toadd_ul);
		toadd_ul.release();
		
		channel_weights(total_channels) = ab_weight / 2;
		homog_weights(total_channels) = ab_homog_weight / 2;
		
		total_channels++;
		
		toadd_ul.reset(new matrix<unsigned long>());
		*toadd_ul = lib_image::quantize_values(lab[2], num_bins);
		channels.add(*toadd_ul);
		toadd_ul.release();
		
		channel_weights(total_channels) = ab_weight / 2;
		homog_weights(total_channels) = ab_homog_weight / 2;
		
		total_channels++;
	}
	
	//texture channel
	if (text_weight < 1e-9 && text_homog_weight < 1e-9) {
		cout << "ignoring texture" << endl;
	} else {
		toadd_ul.reset(new matrix<unsigned long>(t_assign));
		channels.add(*toadd_ul);
		toadd_ul.release();
	
		channel_weights(total_channels) = text_weight;
		homog_weights(total_channels) = text_homog_weight;
		
		total_channels++;
	}
	
	//eigenvector channels
	if (eigvect_weight < 1e-9 && eigvect_homog_weight < 1e-9) {
		cout << "Ignoring eigvects" << endl;
	} else {
		//eigen vector channels, only add if weight is positive
		for (unsigned long i = 0; i < num_eigvects; i++) {
			toadd_ul.reset(new matrix<unsigned long>());
			*toadd_ul = lib_image::quantize_values((*eigvects)[i], num_bins);
			channels.add(*toadd_ul);
			toadd_ul.release();
			
			channel_weights(total_channels) = eigvect_weight / num_eigvects;
			homog_weights(total_channels) = eigvect_homog_weight / num_eigvects;
			
			total_channels++;
		}
	}
	
	//resize parameters to correct size
	channel_weights.resize(total_channels);
	homog_weights.resize(total_channels);
	
	//TODO mirror border
	
	// BUILD OUTPUT PARAMETERS FOR PJ
	//outputs... TODO assume all are same
	matrix<> pjs(im_size_x, im_size_y);
	matrix< array<double> > jangles(im_size_x, im_size_y);
	matrix<> init_ests(im_size_x, im_size_y);
	
	try {
		lib_image::compute_pj_exp(
			channels, poschan, poschan_ori,
			channel_weights, homog_weights, pos_channel_weight,
			threshold_rad_support, pos_channel_threshold,
			n_slices, n_oris_gpb, rad_support, rad_inner_support,
			min_n_angles, max_n_angles, num_bins, norm_power,
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










