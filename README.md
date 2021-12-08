# C++ Library for Gaussian Quadrature with Singular weights

## Installation


To use this header-ony library, clone this repo and run the following commands in the directory you're cloned to:
### Setup with CMake
```
$ cd singular-gauss-quad
$ mkdir build
$ cd build
$ cmake ..
$ sudo cmake --build . --target install
```

The `sudo` is very important, as files will be crated & require admin privileges.

### Setup by including

As this is a header-only library, you can also dowload ``gauss-jacobi.hpp`` and put it into your 'include' file, for exmaple. However, it is important to note the following:
### Dependencies
GQJacobi is a library that relies on Eigen, which you can download by running;
```
brew install eigen
```
If you are building with CMake, this isn't a problem. However, if you are simply including into your `include` directory, you will have to set your compiler path to wherever you have installed Eigen.


## Building

### With CMake
Your master CMake will need to include the following lines

```
find_package(GQJacobi)
target_link_libraries(yourprojectname PUBLIC GQJacobi::GQJacobi)
```

for `QJacobi` to work. The rest is already handled.

### As an include
You will need to configure your project in such a way that the `gauss-jacobi.hpp` file can properly include `Eigen`.

Once `GQJacobi` is installed and incuded, import is by typing:
```
#include "GQJacobi/gauss-jacobi.hpp"
```
