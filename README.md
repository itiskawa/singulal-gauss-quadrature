# SingGQ : a C++ Library for Gaussian Quadrature with Singular weights

### Dependencies
SingGQ is a library that relies on Eigen, which you can download by running;
```
brew install eigen
```
If you are building with CMake, this isn't a problem. However, if you are simply including into your `include` directory, you will have to set your compiler path to wherever you have installed Eigen.

You will also need to have `Boost` installed, which you can do by typing:
```
brew install boost
```


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
If you get the message `Could NOT find Boost (missing: Boost_DIR)` or `Could NOT find Eigen3 (missing: Eigen3_DIR)` when installing, then go back to the previous step and install the dependencies.
The `sudo` is very important, as files will be crated & require admin privileges.

### Setup by including

As this is a header-only library, you can also dowload ``gauss-jacobi.hpp`` and put it into your 'include' file, for exmaple. However, it is important to note the following:



## Building

### With CMake
Your root CMake will need to include the following lines

```
find_package(SingGQ)
target_link_libraries(yourprojectname PUBLIC SingGQ::SingGQ)
```

for `SingGQ` to work. The rest is already handled (dependencies and such).

### As an include
You will need to configure your project in such a way that the `gauss-jacobi.hpp` file can properly include `Eigen`.


**NOTE!** _You may need to add the include path to Eigen in your environment! The code compiles and runs, but the IDE or editor may fuss about, thus making the experience uncomfortable._

## Usage

Once `SingGQ` is installed and incuded, import is by typing:
```
#include "SingGQ/gauss-jacobi.hpp"
```

### Example
Say you want to compute the integral of `cos` with respect to a Chebyshev weight function
````
//instantiate quadrature rule: 10 nodes
GQJacobi::GaussChebyshevRule<double> gch(10, -1); // equivalent to GaussJacobiRule<double> gch(10, -0.5, -0.5)

//compute quadrature
double cosQuadCheb = gch(cos<double>);
````

Lambda notation is also supported

````
double cosQuadCheb = gch([](double x){return cos(x);})
````

Integral computation over different intervals is possible as well (say over ]0,2[)

````
double cosQuadCheb2 = gch(cos<double>, 0,2);
````
