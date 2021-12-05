#include <iostream>
using namespace std;

#include <eigen3/Eigen/Dense>
#include "NumCpp.hpp"
#include <cstdlib>


int main(){

    Eigen::MatrixXd a = Eigen::MatrixXd::Random(3,3);
    cout << a << endl;

    return 0;
}