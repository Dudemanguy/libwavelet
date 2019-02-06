# libwavelet
libwavelet is a very simple library for providing 1 and 2 dimensional wavelet transformation functions. Currently, only the Morlet mother wavelet function is provided. However, other wavelet types can be generated on their own and used with the transform function. Currently both the radix-2 and radix-4 algorithms are implemented. They can be selected with the `wavelet_transform` function.

## Installation
libwavelet uses meson for easy building and installation. After checking out the source, navigate to the directory, and then simply run.
```
mkdir build
meson build
ninja -C build
sudo ninja -C build install
```

## Usage
libwavelet can be incorporated into your code in the usual way (i.e. `#include <wavelet.h>` and `gcc -lwavelet`). The [header](https://github.com/Dudemanguy911/libwavelet/blob/master/include/wavelet.h) shows available functions and brief descriptions. All `fft`, `ifft`, and `wavelet_transform` functions require the length of the array to be exactly equal to a power of two or four (or else an error is shown).

Consider creating a 1 dimensional Morlet transform with a bandwidth of 10 and a central frequency of 5 hz from a signal sampled at 2000 hz. First, create a double complex array, `z` (the imaginary part can be all zeros) that contains your signal. Then create a double array, `time`,  containing the time values of the signal.
```
struct wavelet wave = init_wavelet("morlet", 10, 5, 2000);
double complex *mother = wavelet_mother1(wave, time);
double complex *transform = wavelet_transform1(wave, mother, z, 2);
```
Note that `wavelet_transform1` and `wavelet_transform2` will work with your own custom mother wavelet array. The value of `2` here uses the radix-2 algorithm as part of the transform. To use the radix-4 instead, simply put `4`.

## License
GPLv2 or later.
