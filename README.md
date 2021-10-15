


# Mat3D #
Mat3D software, a toolkit for fast material modeling and design powered by Fast Fourier Transform (FFT) methods.

## Compile Directions ##

### Dependencies ###

Dependencies for the project are:

* FFTW3
* VTK 6.x or 7.x
* Opt++
* NLOPT 

You also need these header libraries:

* Eigen
* CImg


### Compile & Build ###

Make sure dependencies are installed properly. Then to compile in setup mode

```

mkdir build && cd build && cmake ..



compiles the software by

```

make -j

