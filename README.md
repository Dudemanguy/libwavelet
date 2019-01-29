# libwavelet
libwavelet is a very simple library for providing 1 dimensional wavelet transformation functions. Currently, only the morlet wavelet is provided, but other wavelet types can easily be added.

## Installation
libwavelet uses meson for easy building and installation. After checking out the source, navigate to the directory, and then simply run.
```
mkdir build
meson build
ninja -C build
sudo ninja -C build install
```

## Usage
libwavelet can be incorporated into your code in the usual way (i.e. `#include <wavelet.h>` and `gcc -lwavelet`). The [header](https://github.com/Dudemanguy911/libwavelet/blob/master/include/wavelet.h) shows available functions and brief descriptions. Note that if you opt to use `fft` or `ifft` directly, the length of the array must be a power of two (else the computation will never end). The `wavelet_transform` function will pad zeros for you.

Consider creating a Morlet transform with a bandwidth of 10 and a central frequency of 5 hz from a signal sampled at 2000 hz. First, create a `double complex` array, `z` (the imaginary part can be all zeros) that contains your signal.
```
struct wavelet wave = init_wavelet("morlet", 10, 5, 2000);
wavelet_transform(wave, z, k);
```
The transformed wavelet is stored in the `wave.transform` struct member.

## License
GPLv2 or later.
