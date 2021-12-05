# C++ Library for Gaussian Quadrature with Singular weights

## Installation


To use this header-ony library, clone this repo and run the following commands in the directory you're cloned to:
### Setup with CMake
'''
$ cd singular-gauss-quad
$ mkdir build
$ cd build
$ cmake ..
$ sudo cmake --build . --target install
'''

The 'sudo' is very important, as files will be crated & require admin privileges.

### Setup by including
As this is a header-only library, you can also dowload ''gauss-jacobi.hpp'' and put it into your 'include' file, for exmaple. However, it is important to note the following:
### Dependencies
GQJacobi is a library that relies on Eigen, which you can download by running;
'''
brew install eigen
'''
If you are building with CMake, this isn't a problem. However, if you are simply including into your 'include' directory, you will have to set your compiler path to wherever you have installed Eigen.


## Building

To build this 