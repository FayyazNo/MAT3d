
```
___  ___      _   ___________ 
|  \/  |     | | |____ |  _  \
| .  . | __ _| |_    / / | | |
| |\/| |/ _` | __|   \ \ | | |
| |  | | (_| | |_.___/ / |/ / 
\_|  |_/\__,_|\__\____/|___/  

Copyright (c) 2017, SIMSOFTS LLC, All rights reserved.
```

# Mat3D #

This is the serial version of Mat3D software, a toolkit for fast material modeling and design powered by Fast Fourier Transform (FFT) methods.

## Compile Directions ##

### Dependencies ###

Dependencies for the project are:

* FFTW3
* VTK 6.x or 7.x
* Opt++
* NLOPT 


### Compile & Build ###

Make sure dependencies are installed properly. Then to compile in setup mode

```
mkdir build && cd build && cmake ..
```

compiles the software by

```
make -j
```

and then test the build by

```
make test
```

make sure all the tests are passing and then install

```
make install
```

### Recompiling ###

In case you made changes in the source code, recompile the project in the build folder by

```
make
```

### Installing from Package ###

This feature is currently under testing and development. The pacakge is sucessfully installed on Ubuntu 16.04 LTS. To install pacakge:

```
sudo dpkg -i ./MATFFT3D-X.Y.Z-Linudx.deb
```

sucessfully installed software can be permanently removed by

```
sudo dpkg --remove MATFFT3D
```
