#pragma once
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cstddef>
#include <type_traits>
#include <complex>
#include <vector>
#include <cstddef>
#include <type_traits>


/**
 * @class (abstract)
 * @author V.B. (alias @itiskawa)
 * @brief An abstract class that defines common functions to different quadrature solvers
 * It has been chosen to not contain fields nor evaluation operators, due to subtleties in subclasses (notably GaussLogRule<T>)
 * @param T : the data type of the quadrature (float, double etc...)
 * @date February 2022
 */
template<class T>
class GaussRule{

    public : 
    /**
     * @brief Computes the recurrence relation coefficients of the monic Jacobi Polynomials
     * The methods mimics r_jacobi.m of W. Gautschi in his MATLAB Suite for Orthognonal Polynomials
     * available at: https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
     * 
     * @param n : desired number of coefficients, should match the number of quadrature nodes / weights for a Gauss-Jacobi quadrature rule
     * @param a : the a parameter of the Jacobi weight function
     * @param b : the b parameter of the Jacobi weight function
     * @return a nx2 array with the recurrence relation coefficients in the columns
    */ 
    virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> c_jacobi(std::size_t n, double a, double b) {

        // Empty coefficient matrix: alpha_k and beta_k are stored in columns
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> coeffs = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, 2);
        

        // Following Gautschi's method
        double a0 = (b-a)/(a+b+2);
        double b0 = pow(2, (a+b+1))*((tgamma(a+1)*tgamma(b+1))/tgamma(a+b+2));

        // Special first degree case
        if(n == 1) {
            coeffs(0, 0) = a0;
            coeffs(0, 1) = b0;
            return coeffs;
        }

        // Computing / storing the alpha coefficients
        coeffs(0, 0) = a0;
        Eigen::Vector<T,Eigen::Dynamic> deg = Eigen::Vector<T, Eigen::Dynamic>::LinSpaced(n-1,1, n-1);
        Eigen::Vector<T,Eigen::Dynamic> ndeg = (2*deg)+(Eigen::Vector<T, Eigen::Dynamic>::Ones(n-1)*(a+b));

        for(int i = 1; i < n; i++){
            coeffs(i, 0) = (pow(b,2)-pow(a,2))/(ndeg[i-1]*(ndeg[i-1]+2));
        }


        // Computimg / storing the beta coefficients
        coeffs(0, 1) = b0;
        coeffs(1, 1) = 4*(a+1)*(b+1)/(pow(a+b+2, 2)*(a+b+3));
        for(int i = 2; i < n; i++){
            coeffs(i, 1) = (4*(i*(i+b)*(i+a)*(i+a+b))) / (pow(ndeg[i-1], 2)*(ndeg[i-1]-1)*(ndeg[i-1]+1));
        }

        return coeffs;
    }    


    /**
     * @brief rearranges an array of recurrence relation coefficients into a self-adjoint matrix for Golub-Welsch solve
     * The alpha coefficients are placed in the diagonal, and the n-1 first beta coefficients' square roots are placed tridiagonally
     * 
     * example with n=4:
     * 
     *      | a_0 c_1 0    0   0  |
     *      | c_1 a_1 b_2  0   0  |
     * J =  |  0  c_2 a_3 c_2  0  |   where c_i = sqrt(b_i)
     *      |  0   0  c_3 a_4 c_3 |
     * 
     * 
     * @param coeffs: nx2 array of recurrence relation coefficients, n: number of desired coefficients
     * @param n : desired matrix size, should match the number of coefficients
     * @return nxn self-adjoint tridiagonal matrix 
    */
    virtual Eigen::Matrix<T, Eigen::Dynamic,Eigen::Dynamic> tridiagCoeffs(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> coeffs, std::size_t n) {

        // SIZE CHECK
        assert(coeffs.rows() == n);
        assert(coeffs.cols() == 2);
        
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tridiag = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, n);

        // Setting alpha0 in top left corner
        tridiag(0, 0) = coeffs(0, 0); 

        // Setting tridiagonal coefficients
        for(int i = 1; i < n; i++){
            tridiag(i, i) = coeffs(i, 0);
            tridiag(i,i-1) = sqrt(coeffs(i, 1));
            tridiag(i-1, i) = sqrt(coeffs(i, 1));
        }

        return tridiag;
    }

}; // GaussRule<T>


