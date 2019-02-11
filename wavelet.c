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

void fft1_radix2(double complex *z, int len, double complex *tmp) {
	if (len > 1) {
		double complex w, *ze, *zo;
		ze = tmp;
		zo = tmp+len/2;
		for (int k = 0; k < len/2; ++k) {
			ze[k] = z[2*k];
			zo[k] = z[2*k+1];
		}
		fft1_radix2(ze, len/2, z);
		fft1_radix2(zo, len/2, z);
		for (int m = 0; m < len/2; ++m) {
			w = cexp((-2*M_PI*I*m)/(double)len);
			z[m] = ze[m] + w*zo[m];
			z[m+len/2] = ze[m] - w*zo[m];
		}
	}
}

void fft1_radix4(double complex *z, int len, double complex *tmp) {
	if (len > 1) {
		double complex w1, w2, w3, *z0, *z1, *z2, *z3;
		z0 = tmp;
		z1 = tmp+len/4;
		z2 = tmp+len/2;
		z3 = tmp+(3*len)/4;
		for (int k = 0; k < len/4; ++k) {
			z0[k] = z[4*k];
			z1[k] = z[4*k+1];
			z2[k] = z[4*k+2];
			z3[k] = z[4*k+3];
		}
		fft1_radix4(z0, len/4, z);
		fft1_radix4(z1, len/4, z);
		fft1_radix4(z2, len/4, z);
		fft1_radix4(z3, len/4, z);
		for (int m = 0; m < len/4; ++m) {
			w1 = cexp((-2*M_PI*I*(m))/(double)len);
			w2 = cexp((-2*M_PI*I*(2*m))/(double)len);
			w3 = cexp((-2*M_PI*I*(3*m))/(double)len);
			z[m] = z0[m] + w1*z1[m] + w2*z2[m] + w3*z3[m];
			z[m+len/4] = z0[m] - I*w1*z1[m] - w2*z2[m] + I*w3*z3[m];
			z[m+len/2] = z0[m] - w1*z1[m] + w2*z2[m] - w3*z3[m];
			z[m+(3*len)/4] = z0[m] + I*w1*z1[m] - w2*z2[m] - I*w3*z3[m];
		}
	}
}

void fft2_radix2(double complex **z, int len, double complex **tmp, int width) {
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
		fft2_radix2(ze, len/2, z, width);
		fft2_radix2(zo, len/2, z, width);
		for (int m = 0; m < len/2; ++m) {
			w = cexp((-2*M_PI*I*m)/(double)len);
			for (int n = 0; n < width; ++n) {
				z[m][n] = ze[m][n] + w*zo[m][n];
				z[m+len/2][n] = ze[m][n] - w*zo[m][n];
			}
		}
	}
}

void fft2_radix4(double complex **z, int len, double complex **tmp, int width) {
	if (len > 1) {
		double complex w1, w2, w3, **z0, **z1, **z2, **z3;
		z0 = tmp;
		z1 = tmp+len/4;
		z2 = tmp+len/2;
		z3 = tmp+(3*len)/4;
		for (int k = 0; k < len/4; ++k) {
			for (int n = 0; n < width; ++n) {
				z0[k][n] = z[4*k][n];
				z1[k][n] = z[4*k+1][n];
				z2[k][n] = z[4*k+2][n];
				z3[k][n] = z[4*k+3][n];
			}
		}
		fft2_radix4(z0, len/4, z, width);
		fft2_radix4(z1, len/4, z, width);
		fft2_radix4(z2, len/4, z, width);
		fft2_radix4(z3, len/4, z, width);
		for (int m = 0; m < len/4; ++m) {
			w1 = cexp((-2*M_PI*I*(m))/(double)len);
			w2 = cexp((-2*M_PI*I*(2*m))/(double)len);
			w3 = cexp((-2*M_PI*I*(3*m))/(double)len);
			for (int n = 0; n < width; ++n) {
				z[m][n] = z0[m][n] + w1*z1[m][n] + w2*z2[m][n] + w3*z3[m][n];
				z[m+len/4][n] = z0[m][n] - I*w1*z1[m][n] - w2*z2[m][n] + I*w3*z3[m][n];
				z[m+len/2][n] = z0[m][n] - w1*z1[m][n] + w2*z2[m][n] - w3*z3[m][n];
				z[m+(3*len)/4][n] = z0[m][n] + I*w1*z1[m][n] - w2*z2[m][n] - I*w3*z3[m][n];
			}
		}
	}
}

void fft1(double complex *z, int len, double complex *tmp, int radix) {
	if (radix == 2) {
		fft1_radix2(z, len, tmp);
	}
	if (radix == 4) {
		fft1_radix4(z, len, tmp);
	}
}

void fft2(double complex **z, int len, double complex **tmp, int width, int radix) {
	if (radix == 2) {
		fft2_radix2(z, len, tmp, width);
	}
	if (radix == 4) {
		fft2_radix4(z, len, tmp, width);
	}
}

void ifft1_rrad2(double complex *z, int len, double complex *tmp) {
	if (len > 1) {
		int k, m;
		double complex w, *zo, *ze;
		ze = tmp;
		zo = tmp+len/2;
		for (k = 0; k < len/2; ++k) {
			ze[k] = z[2*k];
			zo[k] = z[2*k+1];
		}
		ifft1_rrad2(ze, len/2, z);
		ifft1_rrad2(zo, len/2, z);
		for (m = 0; m < len/2; ++m) {
			w = cexp((2*M_PI*I*m)/(double)len);
			z[m] = ze[m] + w*zo[m];
			z[m+len/2] = ze[m] - w*zo[m];
		}
	}
}

void ifft1_rrad4(double complex *z, int len, double complex *tmp) {
	if (len > 1) {
		double complex w1, w2, w3, *z0, *z1, *z2, *z3;
		z0 = tmp;
		z1 = tmp+len/4;
		z2 = tmp+len/2;
		z3 = tmp+(3*len)/4;
		for (int k = 0; k < len/4; ++k) {
			z0[k] = z[4*k];
			z1[k] = z[4*k+1];
			z2[k] = z[4*k+2];
			z3[k] = z[4*k+3];
		}
		ifft1_rrad4(z0, len/4, z);
		ifft1_rrad4(z1, len/4, z);
		ifft1_rrad4(z2, len/4, z);
		ifft1_rrad4(z3, len/4, z);
		for (int m = 0; m < len/4; ++m) {
			w1 = cexp((2*M_PI*I*(m))/(double)len);
			w2 = cexp((2*M_PI*I*(2*m))/(double)len);
			w3 = cexp((2*M_PI*I*(3*m))/(double)len);
			z[m] = z0[m] + w1*z1[m] + w2*z2[m] + w3*z3[m];
			z[m+len/4] = z0[m] + I*w1*z1[m] - w2*z2[m] - I*w3*z3[m];
			z[m+len/2] = z0[m] - w1*z1[m] + w2*z2[m] - w3*z3[m];
			z[m+(3*len)/4] = z0[m] - I*w1*z1[m] - w2*z2[m] + I*w3*z3[m];
		}
	}
}

void ifft2_rrad2(double complex **z, int len, double complex **tmp, int width) {
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
		ifft2_rrad2(ze, len/2, z, width);
		ifft2_rrad2(zo, len/2, z, width);
		for (int m = 0; m < len/2; ++m) {
			w = cexp((2*M_PI*I*m)/(double)len);
			for (int n = 0; n < width; ++n) {
				z[m][n] = ze[m][n] + w*zo[m][n];
				z[m+len/2][n] = ze[m][n] - w*zo[m][n];
			}
		}
	}
}

void ifft2_rrad4(double complex **z, int len, double complex **tmp, int width) {
	if (len > 1) {
		double complex w1, w2, w3, **z0, **z1, **z2, **z3;
		z0 = tmp;
		z1 = tmp+len/4;
		z2 = tmp+len/2;
		z3 = tmp+(3*len)/4;
		for (int k = 0; k < len/4; ++k) {
			for (int n = 0; n < width; ++n) {
				z0[k][n] = z[4*k][n];
				z1[k][n] = z[4*k+1][n];
				z2[k][n] = z[4*k+2][n];
				z3[k][n] = z[4*k+3][n];
			}
		}
		ifft2_rrad4(z0, len/4, z, width);
		ifft2_rrad4(z1, len/4, z, width);
		ifft2_rrad4(z2, len/4, z, width);
		ifft2_rrad4(z3, len/4, z, width);
		for (int m = 0; m < len/4; ++m) {
			w1 = cexp((2*M_PI*I*(m))/(double)len);
			w2 = cexp((2*M_PI*I*(2*m))/(double)len);
			w3 = cexp((2*M_PI*I*(3*m))/(double)len);
			for (int n = 0; n < width; ++n) {
				z[m][n] = z0[m][n] + w1*z1[m][n] + w2*z2[m][n] + w3*z3[m][n];
				z[m+len/4][n] = z0[m][n] + I*w1*z1[m][n] - w2*z2[m][n] - I*w3*z3[m][n];
				z[m+len/2][n] = z0[m][n] - w1*z1[m][n] + w2*z2[m][n] - w3*z3[m][n];
				z[m+(3*len)/4][n] = z0[m][n] - I*w1*z1[m][n] - w2*z2[m][n] + I*w3*z3[m][n];
			}
		}
	}
}

void ifft1_radix2(double complex *z, int len, double complex *tmp) {
	ifft1_rrad2(z, len, tmp);
	for (int i = 0; i < len; ++i) {
		z[i] = (1/(double)len)*z[i];
	}
}

void ifft1_radix4(double complex *z, int len, double complex *tmp) {
	ifft1_rrad4(z, len, tmp);
	for (int i = 0; i < len; ++i) {
		z[i] = (1/(double)len)*z[i];
	}
}

void ifft1(double complex *z, int len, double complex *tmp, int radix) {
	if (radix == 2) {
		ifft1_radix2(z, len, tmp);
	}
	if (radix == 4) {
		ifft1_radix4(z, len, tmp);
	}
}

void ifft2_radix2(double complex **z, int len, double complex **tmp, int width) {
	ifft2_rrad2(z, len, tmp, width);
	for (int i = 0; i < len; ++i) {
		for (int j = 0; j < width; ++j) {
			z[i][j] = (1/(double)len)*z[i][j];
		}
	}
}

void ifft2_radix4(double complex **z, int len, double complex **tmp, int width) {
	ifft2_rrad4(z, len, tmp, width);
	for (int i = 0; i < len; ++i) {
		for (int j = 0; j < width; ++j) {
			z[i][j] = (1/(double)len)*z[i][j];
		}
	}
}

void ifft2(double complex **z, int len, double complex **tmp, int width, int radix) {
	if (radix == 2) {
		ifft2_radix2(z, len, tmp, width);
	}
	if (radix == 4) {
		ifft2_radix4(z, len, tmp, width);
	}
}

int power_of_two_or_four(int len) {
	double x = log2(len);
	double y = log(len)/log(4);
	if (x == (int)x) {
		return 1;
	} else  if (y == (int)y) {
		return 1;
	} else {
		return 0;
	}
}

struct wavelet init_wavelet(char *type, double bwidth, double cfq, double srate, int len, int width) {
	struct wavelet wave;
	if (!len) {
		printf("Error, wavelet object has length 0. Please reinitialize the"
				"wavelet with a nonzero power of two or four.\n");
		wave.error = 1;
		return wave;
	}
	if (!power_of_two_or_four(len)) {
		printf("Error, wavelet object must have a length equal to a nonzero"
				"power of two or four.\n");
		wave.error = 1;
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

double complex *wavelet_transform1(struct wavelet *wave, double complex *mother, double complex *z, int radix) {
	if (radix != 2 || radix != 4) {
		radix = 2;
	}
	double complex *signal = (double complex *)malloc(sizeof(double complex)*wave->len);
	memcpy(signal, z, wave->len*sizeof(&z));
	double complex *tmp = (double complex *)malloc(sizeof(double complex)*wave->len);
	fft1(signal, wave->len, tmp, radix);
	double complex *mother_tmp = (double complex *)malloc(sizeof(double complex)*wave->len);
	memcpy(mother_tmp, mother, wave->len*sizeof(&mother));
	fft1(mother_tmp, wave->len, tmp, radix);
	double complex *transform = (double complex *)malloc(sizeof(double complex)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		transform[i] = complex_multiply(signal[i], mother_tmp[i]);
	}
	ifft1(transform, wave->len, tmp, radix);
	free(mother_tmp);
	free(signal);
	free(tmp);
	return transform;
}

double complex **wavelet_transform2(struct wavelet *wave, double complex **mother, double complex **z, int radix) {
	if (radix != 2 || radix != 4) {
		radix = 2;
	}
	double complex **signal = (double complex **)malloc(sizeof(double complex *)*wave->len);
	memcpy(signal, z, wave->len*sizeof(&z));
	double complex **tmp = (double complex **)malloc(sizeof(double complex *)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		tmp[i] = (double complex *)malloc(sizeof(double complex)*wave->width);
	}
	fft2(signal, wave->len, tmp, wave->width, radix);
	double complex **mother_tmp = (double complex **)malloc(sizeof(double complex *)*wave->len);
	memcpy(mother_tmp, mother, wave->len*sizeof(&mother));
	fft2(mother_tmp, wave->len, tmp, wave->width, radix);
	double complex **transform = (double complex **)malloc(sizeof(double complex *)*wave->len);
	for (int i = 0; i < wave->len; ++i) {
		transform[i] = (double complex *)malloc(sizeof(double complex)*wave->width);
	}
	for (int i = 0; i < wave->width; ++i) {
		for (int j = 0; j < wave->len; ++j) {
			transform[j][i] = complex_multiply(signal[j][i], mother_tmp[j][i]);
		}
	}
	ifft2(transform, wave->len, tmp, wave->width, radix);
	for (int i = 0; i < wave->len; ++i) {
		free(tmp[i]);
	}
	free(mother_tmp);
	free(signal);
	free(tmp);
	return transform;
}
