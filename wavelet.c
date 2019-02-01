#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wavelet.h"

double complex complex_multiply(double complex v, double complex w) {
	double a = creal(v);
	double b = cimag(v);
	double c = creal(w);
	double d = cimag(w);
	double complex z = (a*c - b*d) + (a*d+b*c)*I;
	return z;
}

void fft1(double complex *z, int len, double complex *tmp) {
	if (len > 1) {
		double complex w, *zo, *ze;
		ze = tmp;
		zo = tmp+len/2;
		for (int k = 0; k < len/2; ++k) {
			ze[k] = z[2*k];
			zo[k] = z[2*k+1];
		}
		fft1(ze, len/2, z);
		fft1(zo, len/2, z);
		for (int m = 0; m < len/2; ++m) {
			w = cexp((-2*M_PI*I*m)/(double)len);
			z[m] = ze[m] + w*zo[m];
			z[m+len/2] = ze[m] - w*zo[m];
		}
	}
}

void fft2(double complex **z, int len, double complex **tmp, int width) {
	if (len > 1) {
		double complex w, **zo, **ze;
		ze = tmp;
		zo = tmp+len/2;
		for (int k = 0; k < len/2; ++k) {
			for (int n = 0; n < width; ++n) {
				ze[k][n] = z[2*k][n];
				zo[k][n] = z[2*k+1][n];
			}
		}
		fft2(ze, len/2, z, width);
		fft2(zo, len/2, z, width);
		for (int m = 0; m < len/2; ++m) {
			w = cexp((-2*M_PI*I*m)/(double)len);
			for (int n = 0; n < width; ++n) {
				z[m][n] = ze[m][n] + w*zo[m][n];
				z[m+len/2][n] = ze[m][n] - w*zo[m][n];
			}
		}
	}
}

void ifft1_rec(double complex *z, int len, double complex *tmp) {
	if (len > 1) {
		int k, m;
		double complex w, *zo, *ze;
		ze = tmp;
		zo = tmp+len/2;
		for (k = 0; k < len/2; ++k) {
			ze[k] = z[2*k];
			zo[k] = z[2*k+1];
		}
		ifft1_rec(ze, len/2, z);
		ifft1_rec(zo, len/2, z);
		for (m = 0; m < len/2; ++m) {
			w = cexp((2*M_PI*I*m)/(double)len);
			z[m] = ze[m] + w*zo[m];
			z[m+len/2] = ze[m] - w*zo[m];
		}
	}
}

void ifft2_rec(double complex **z, int len, double complex **tmp, int width) {
	if (len > 1) {
		double complex w, **zo, **ze;
		ze = tmp;
		zo = tmp+len/2;
		for (int k = 0; k < len/2; ++k) {
			for (int n = 0; n < width; ++n) {
				ze[k][n] = z[2*k][n];
				zo[k][n] = z[2*k+1][n];
			}
		}
		ifft2_rec(ze, len/2, z, width);
		ifft2_rec(zo, len/2, z, width);
		for (int m = 0; m < len/2; ++m) {
			w = cexp((2*M_PI*I*m)/(double)len);
			for  (int n = 0; n < width; ++n) {
				z[m][n] = ze[m][n] + w*zo[m][n];
				z[m+len/2][n] = ze[m][n] - w*zo[m][n];
			}
		}
	}
}

void ifft1(double complex *z, int len, double complex *tmp) {
	ifft1_rec(z, len, tmp);
	for (int i = 0; i < len; ++i) {
		z[i] = (1/(double)len)*z[i];
	}
}

void ifft2(double complex **z, int len, double complex **tmp, int width) {
	ifft2_rec(z, len, tmp, width);
	for (int i = 0; i < len; ++i) {
		for (int j = 0; j < width; ++j) {
			z[i][j] = (1/(double)len)*z[i][j];
		}
	}
}

int power_of_two(int len) {
	double x = log2(len);
	if (x == (int)x) {
		return 1;
	} else {
		return 0;
	}
}

struct wavelet init_wavelet(char *type, double bwidth, double cfq, double srate, int len, int width) {
	struct wavelet wave;
	if (!len) {
		printf("Error, wavelet object has length 0. Please reinitialize the"
				"wavelet with a nonzero power of two.\n");
		return wave;
	}
	if (!power_of_two(len)) {
		printf("Error, wavelet object must have a length equal to a nonzero"
				"power of two.\n");
		return wave;
	}
	wave.bwidth = bwidth;
	wave.cfq = cfq;
	wave.srate = srate;
	wave.len = len;
	wave.width = width;
	if (strcmp(type, "morlet") == 0) {
		wave.type = MORL;
	}
	return wave;
}

double complex *wavelet_mother1(struct wavelet *wave, double *time) {
	if (!wave) {
		printf("Error, wavelet object must be initialized first with init_wavelet.\n");
		return 0;
	}
	double complex *mother = (double complex *)malloc(sizeof(double complex)*wave->len);
	if (wave->type == MORL) {
		for (int i = 0; i < wave->len; ++i) {
			mother[i] = (1/sqrt(M_PI*wave->bwidth))*cexp((-1*pow(time[i],2))/wave->bwidth)*
				cexp(-2*M_PI*I*wave->cfq*time[i]);
		}
	}
	return mother;
}

double complex **wavelet_mother2(struct wavelet *wave, double *time) {
	if (!wave) {
		printf("Error, wavelet object must be initialized first with init_wavelet.\n");
		return 0;
	}
	double complex **mother = (double complex **)malloc(sizeof(double complex)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		mother[i] = (double complex *)malloc(sizeof(double complex)*wave->width);
	}
	if (wave->type == MORL) {
		for (int i = 0; i < wave->width; ++i) {
			for (int j = 0; j < wave->len; ++j) {
				if ((time[j] == 0) && (j > 0)) {
					mother[j][i] = 0 + 0*I;
				} else {
					mother[j][i] = (1/sqrt(M_PI*wave->bwidth))*cexp((-1*pow(time[j],2))/wave->bwidth)*
						cexp(-2*M_PI*I*wave->cfq*time[j]);
				}
			}
		}
	}
	return mother;
}

double complex *wavelet_transform1(struct wavelet *wave, double complex *mother, double complex *z) {
	//save a copy of the input
	double complex *signal = (double complex *)malloc(sizeof(double complex)*wave->len);
	memcpy(signal, z, wave->len*sizeof(&z));
	double complex *tmp = (double complex *)malloc(sizeof(double complex)*wave->len);
	fft1(signal, wave->len, tmp);
	//save a copy of the mother wavelet
	double complex *mother_tmp = (double complex *)malloc(sizeof(double complex)*wave->len);
	memcpy(mother_tmp, mother, wave->len*sizeof(&mother));
	//transform mother wavelet to frequency space and multiply element-wise
	fft1(mother_tmp, wave->len, tmp);
	double complex *transform = (double complex *)malloc(sizeof(double complex)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		transform[i] = complex_multiply(signal[i], mother_tmp[i]);
	}
	//transform back to time space
	ifft1(transform, wave->len, tmp);
	free(signal);
	free(tmp);
	free(mother_tmp);
	return transform;
}

double complex **wavelet_transform2(struct wavelet *wave, double complex **mother, double complex **z) {
	//save a copy of the input
	double complex **signal = (double complex **)malloc(sizeof(double complex *)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		signal[i] = (double complex *)malloc(sizeof(double complex)*wave->width);
	}
	memcpy(signal, z, wave->len*sizeof(&z));
	double complex **tmp = (double complex **)malloc(sizeof(double complex *)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		tmp[i] = (double complex *)malloc(sizeof(double complex)*wave->width);
	}
	fft2(signal, wave->len, tmp, wave->width);
	//save a copy of the mother wavelet
	double complex **mother_tmp = (double complex **)malloc(sizeof(double complex *)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		mother_tmp[i] = (double complex *)malloc(sizeof(double complex)*wave->width);
	}
	memcpy(mother_tmp, mother, wave->len*sizeof(&mother));
	//transform mother wavelet to frequency space and multiple element-wise
	fft2(mother_tmp, wave->len, tmp, wave->width);
	double complex **transform = (double complex **)malloc(sizeof(double complex *)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		transform[i] = (double complex *)malloc(sizeof(double complex)*wave->width);
	}
	for (int i = 0; i < wave->width; ++i) {
		for (int j = 0; j < wave->len; ++j) {
			transform[j][i] = complex_multiply(signal[j][i], mother_tmp[j][i]);
		}
	}
	//transform back to time space
	ifft2(transform, wave->len, tmp, wave->width);
	for (int i = 0; i < wave->len; ++i) {
		free(signal[i]);
		free(mother_tmp[i]);
		free(tmp[i]);
	}
	free(signal);
	free(mother_tmp);
	free(tmp);
	return transform;
}
