/**
 * Executes the fast fourier transform on the double complex array, z. 
 * z must be padded with zeros to be exactly the length of a power of 2. n is
 * merely the length of v. tmp is a temporary array with allocated memory of size
 * n for computational purposes. fft1 is the one dimensional transform and fft2
 * is the two dimensional transform. For fft2, an extra parameter specifying
 * the width of the array is required.
 */
void fft1(double complex *z, int len, double complex *tmp);
void fft2(double complex **z, int len, double complex **tmp, int width);

/**
 * Executes the inverse fast fourier transform on the double complex array, z. 
 * z must be padded with zeros to be exactly the length of a power of 2. n is
 * merely the length of v. tmp is a temporary array with allocated memory of size
 * n for computational purposes. ifft1 is the one dimensional transform and ifft2
 * is the two dimensional transform. For ifft2, an extra parameter specifying
 * the width of the array is required.
 */
void ifft1(double complex *z, int len, double complex *tmp);
void ifft2(double complex **z, int len, double complex **tmp, int width);

/**
 * Initiates the wavelet object. Currently, the only supported type is "morlet", 
 * but more may be added in the future. bdwidth is the bandwith parameter of a
 * wavelet which dilates/compresses the signal. cfq is the central frequency of 
 * the wavelet. srate is the sample rate of the input signal. len is required
 * and must be equal to a power of two or else the function will return an error.
 * For a one dimensional transform, simply pass width as 0. For two dimensions, 
 * make width equal to the width of your input signal.
 */
struct wavelet init_wavelet(char *type, double bwidth, double cfq, double srate, int len, int width);

/**
 * Returns either a one or two dimensional array containing entries of the mother 
 * wavelet set by init_wavelet. Currently only "morlet" is supported. The 
 * wavelet_mother functions also require an array of doubles containing the time 
 * axis that corresponds to the signal. Note that the time array must be equal 
 * to wave->len.
 */
double complex *wavelet_mother1(struct wavelet *wave, double *time);
double complex **wavelet_mother2(struct wavelet *wave, double *time);

/**
 * Returns either the one or two dimensional wavelet transformation.
 * The wavelet object must be initialized first. A morlet mother wavelet may be 
 * obtained via the wavelet_mother functions. However, using a custom mother
 * arrray generated on your own is also valid. 
 */
double complex *wavelet_transform1(struct wavelet *wave, double complex *mother, double complex *z);
double complex **wavelet_transform2(struct wavelet *wave, double complex **mother, double complex **z);

enum wavelet_type {
	MORL,
};

struct wavelet {
	enum wavelet_type *type;
	double bwidth;
	double cfq;
	double srate;
	int len;
	int width;
};
