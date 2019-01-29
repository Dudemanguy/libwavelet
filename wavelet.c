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

void fft(double complex *z, int n, double complex *tmp) {
	if (n > 1) {
		int k, m;
		double complex w, *zo, *ze;
		ze = tmp;
		zo = tmp+n/2;
		for (k = 0; k < n/2; ++k) {
			ze[k] = z[2*k];
			zo[k] = z[2*k+1];
		}
		fft(ze, n/2, z);
		fft(zo, n/2, z);
		for (m = 0; m < n/2; ++m) {
			w = cexp((-2*M_PI*I*m)/(double)n);
			z[m] = ze[m] + w*zo[m];
			z[m+n/2] = ze[m] - w*zo[m];
		}
	}
}

void ifft_rec(double complex *z, int n, double complex *tmp) {
	if (n > 1) {
		int k, m;
		double complex w, *zo, *ze;
		ze = tmp;
		zo = tmp+n/2;
		for (k = 0; k < n/2; ++k) {
			ze[k] = z[2*k];
			zo[k] = z[2*k+1];
		}
		ifft_rec(ze, n/2, z);
		ifft_rec(zo, n/2, z);
		for (m = 0; m < n/2; ++m) {
			w = cexp((2*M_PI*I*m)/(double)n);
			z[m] = ze[m] + w*zo[m];
			z[m+n/2] = ze[m] - w*zo[m];
		}
	}
}

void ifft(double complex *z, int n, double complex *tmp) {
	ifft_rec(z, n, tmp);
	for (int i = 0; i < n; ++i) {
		z[i] = (1/(double)n)*z[i];
	}
}

struct wavelet init_wavelet(char *type, double bwidth, double cfq, double srate) {
	struct wavelet wave;
	wave.bwidth = bwidth;
	wave.cfq = cfq;
	wave.srate = srate;
	if (strcmp(type, "morlet") == 0) {
		wave.type = MORL;
	}
	return wave;
}

void wavelet_transform(struct wavelet wave, double complex *z, int len) {
	int expt = 0;
	while (1) {
		int total = pow(2, expt);
		if (total > len) {
			break;
		} else {
			++expt;
		}
	}
	int n = pow(2, expt);
	double complex *signal = (double complex *)malloc(sizeof(double complex)*n);
	for (int i = 0; i < n; ++i) {
		if (i > len) {
			signal[i] = 0 + 0*I;
		} else {
			signal[i] = z[i];
		}
	}
	//transform signal to frequency space
	double complex *tmp = (double complex *)malloc(sizeof(double complex)*n);
	fft(signal, n, tmp);
	double period = 1/wave.srate;
	double *time = (double *)malloc(sizeof(double)*n);
	for (int i = 0; i < n; ++i) {
		if (i == 0) {
			time[i] = 0;
		} else {
			time[i] = time[i-1] + period;
		}
	}
	//create mother wavelet
	wave.mother = (double complex *)malloc(sizeof(double complex)*n);
	if (wave.type == MORL) {
		for (int i = 0; i < n; ++i) {
			if (i > len) {
				wave.mother[i] = 0;
			} else {
				wave.mother[i] = (1/sqrt(M_PI*wave.bwidth))*cexp((-1*pow(time[i],2))/wave.bwidth)*cexp(-2*M_PI*I*wave.cfq*time[i]);
			}
		}
	}
	double complex *wave_tmp = (double complex *)malloc(sizeof(double complex)*n);
	for (int i = 0; i < n; ++i) {
		wave_tmp[i] = wave.mother[i];
	}
	//transform mother wavelet to frequency space and multiply element-wise
	fft(wave_tmp, n, tmp);
	wave.transform = (double complex *)malloc(sizeof(double complex)*n);
	for (int i = 0; i < n; ++i) {
		wave.transform[i] = complex_multiply(signal[i], wave_tmp[i]);
	}
	//transform back to time space
	ifft(wave.transform, n, tmp);
}
