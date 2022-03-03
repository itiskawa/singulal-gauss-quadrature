#include <iostream>
#include <SingGQ/gauss-jacobi.hpp>
#include <SingGQ/gauss-log.hpp>


using namespace std;



int main(){

    GQJacobi::GaussJacobiRule<double> gq(10, -0.3, -0.7); // 10-point Gauss-Jacobi quadrature rule with alpha = -0.3 and beta = -0.7

    double jacobiSinInteg = gq([](double x){return sin(x);}); // evaluating the weighted sin(x) over ]-1,1[

    double shiftedJacobiSin = gq([](double x){return sin(x);}, 0, 2); // same integral but over ]0,2[

    return 0;
}