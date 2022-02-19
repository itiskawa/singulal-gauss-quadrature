#pragma once

#include <Eigen/Dense>
#include <vector>

namespace GQJacobi{

    struct GaussJacobiRule{
        std::vector<double> nodes;
        std::vector<double> weights;
    };
    
    double gamma_zero(double a, double b) {
        return (pow(2, a+b+1)*tgamma(a+1)*tgamma(b+1))/(tgamma(a+b+2));
    }



    Eigen::MatrixXd c_jacobi(int n, double a, double b) {
        assert(n>1);
        assert(a > -1);
        assert(b > -1);

        // coefficient matrix: alpha and beta stored in columns, goes from 0 to n
        Eigen::MatrixXd coeffs = Eigen::MatrixXd::Zero(n, 2);
        

        // this method follows Gautschi's r_jacobi function
        double a0 = (b-a)/(a+b+2);
        double b0 = pow(2, (a+b+1))*((tgamma(a+1)*tgamma(b+1))/tgamma(a+b+2));

        // special first degree case
        if(n == 1) {
            coeffs(0, 0) = a0;
            coeffs(0, 1) = b0;
            return coeffs;
        }

        // alpha coefficients
        coeffs(0, 0) = a0;
        Eigen::VectorXd deg = Eigen::VectorXd::LinSpaced(n-1,1, n-1);
        Eigen::VectorXd ndeg = (2*deg)+(Eigen::VectorXd::Ones(n-1)*(a+b));

        for(int i = 1; i < n; i++){
            coeffs(i, 0) = (pow(b,2)-pow(a,2))/(ndeg[i-1]*(ndeg[i-1]+2));
        }


        // beta coefficients
        coeffs(0, 1) = b0;
        coeffs(1, 1) = 4*(a+1)*(b+1)/(pow(a+b+2, 2)*(a+b+3));
        for(int i = 2; i < n; i++){
            coeffs(i, 1) = (4*(i*(i+b)*(i+a)*(i+a+b))) / (pow(ndeg[i-1], 2)*(ndeg[i-1]-1)*(ndeg[i-1]+1));
        }
        return coeffs;
    }    



    Eigen::MatrixXd tridiagCoeffs(Eigen::MatrixXd coeffs, int n) {
        // argument is a nx2 matrix
        // SIZE CHECK
        assert(coeffs.rows() == n);
        assert(coeffs.cols() == 2);
        
        Eigen::MatrixXd tridiag = Eigen::MatrixXd::Zero(n, n);

        // setting alpha0 in top left corner
        tridiag(0, 0) = coeffs(0, 0); 

        // setting tridiagonal coefficients
        for(int i = 1; i < n; i++){
            tridiag(i, i) = coeffs(i, 0);
            tridiag(i,i-1) = sqrt(coeffs(i, 1));
            tridiag(i-1, i) = sqrt(coeffs(i, 1));
        }

        return tridiag;
    }

    GaussJacobiRule jacobi_nw(int n, double a, double b) {

        double gamma_0 = gamma_zero(a, b);

        Eigen::MatrixXd coeffs = c_jacobi(n, a, b);
        Eigen::MatrixXd J_n = tridiagCoeffs(coeffs, n);
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solve(J_n); // yields much faster computations of high n
        

        Eigen::VectorXd nodes_v = solve.eigenvalues().real();

        Eigen::MatrixXd eigenvecs = solve.eigenvectors().real();

        Eigen::VectorXd weights_v = Eigen::VectorXd::Zero(n);

        for(int i = 0; i < n; i++){
            weights_v[i] = gamma_0*pow(eigenvecs.col(i).normalized()[0], 2);
        }


        std::vector<double> nodes;
        std::vector<double> weights;


        for(int i = 0; i < n; i++){
            nodes.push_back(nodes_v[i]);
            weights.push_back(weights_v[i]);
        }

        GaussJacobiRule nw = {nodes, weights};
        
        

        return nw;
        
    }



} // namespace GQJacobi