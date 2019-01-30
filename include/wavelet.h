/**
 * Executes the fast fourier transform on the double complex array, z. 
 * z must be padded with zeros to be exactly the length of a power of 2. n is
 * merely the length of v. tmp is a temporary array with allocated memory of size
 * n for computational purposes.
 */
void fft(double complex *z, int n, double complex *tmp);

/**
 * Executes the inverse fast fourier transform on the double complex array, z. 
 * z must be padded with zeros to be exactly the length of a power of 2. n is
 * merely the length of v. tmp is a temporary array with allocated memory of size
 * n for computational purposes.
 */
void ifft(double complex *z, int n, double complex *tmp);

/**
 * Initiates the wavelet object. Currently, the only supported type is "morlet", 
 * but more may be added in the future. bdwidth is the bandwith parameter of a
 * wavelet which dilates/compresses the signal. cfq is the central frequency of 
 * the wavelet. srate is the sample rate of the input signal.
 */
struct wavelet init_wavelet(char *type, double bwidth, double cfq, double srate);

/**
 * Performs the wavelet transform. The wavelet object must be initialized first.
 * wave is, of course, the initialized wavlet object to be used for the transform.
 * z is the double complex array of the input signal and len is the length of z.
 * Note that wavelet_transform handles zero padding for you, so it is not required
 * for z to be a power of 2. The output of the mother wavelet is stored in wave.mother
 * while the output of the transform is stored in wave.transform.
 */
void wavelet_transform(struct wavelet *wave, double complex *z, int len);

/**
 * Destroys memory allocated for the wave object.
 */
void wavelet_destroy(struct wavelet *wave);

enum wavelet_type {
	MORL,
};

struct wavelet {
	enum wavelet_type *type;
	double bwidth;
	double cfq;
	double srate;
	int n;
	double complex *mother;
	double complex *transform;
};
