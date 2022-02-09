#pragma once
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cstddef>
#include <type_traits>
#include <boost/math/special_functions/digamma.hpp>
using namespace Eigen;

namespace GQLog {

    template<class T>
    class GaussLogRule{


        public:

        std::vector<T> nodes;
        std::vector<T> weights;
        std::size_t degree;

        GaussLogRule() = default;

        GaussLogRule(std::size_t n){
            assert(n>1);
            degree = n;

        }

        /*
        * @method 
        * @brief computes the recurrence relation coefficents (alpha_n, beta_n) of the
        * monic polynomials associated to the Jacobi weight function
        */
        Matrix<T, Dynamic, Dynamic> c_jacobi(std::size_t n, double a, double b) {

            // coefficient matrix: alpha and beta stored in columns, goes from 0 to n
            Matrix<T, Dynamic, Dynamic> coeffs = Matrix<T, Dynamic, Dynamic>::Zero(n, 2);
            

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
            VectorXd deg = Vector<T, Dynamic>::LinSpaced(n-1,1, n-1);
            VectorXd ndeg = (2*deg)+(Vector<T, Dynamic>::Ones(n-1)*(a+b));

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


        /*
        * @brief computes the shifted recurrence relation terms
        *
        */
        Matrix<T, Dynamic, Dynamic> shifted_c_log(std::size_t n){
            Matrix<T, Dynamic, Dynamic> abj = Matrix<T, Dynamic, Dynamic>::Zero(n,2);
            Matrix<T, Dynamic, Dynamic> ab = c_jacobi(n, 0, 0);


            //abj(0,0) = (1+ab(0,0))/2.;
            //abj(0,1) = (ab(0,1))/2;

            for(int i = 0; i < n; i++){
                abj(i,0) = (1+ab(i,0))/2.;
            }
            abj(0, 1) = ab(0,1)/2.;
            for(int i = 1; i < n; i++){
                abj(i,1) = ab(i,1)/4.;
            }
            return abj;
        }


        /*
        * @brief computes modified moments
        *
        */
        /*Vector<T, Dynamic> mmom(std::size_t n){
            Vector<T, Dynamic> mom = Vector<T, Dynamic>::Zero(n);

            for(int i = 0; i < n; i++){
                //mom[i] = pow(-1, i)/((i+1)*i);
                mom[i] = 1/(i+1);
            }
            return mom;
        }


        Vector<T, Dynamic> mom_jaclog(std::size_t n){
            Vector<T, Dynamic> mom = Vector<T, Dynamic>::Zero(2*n);
            for(int i = 0; i < 2*n; i++){
                int k = i+1;
                mom[i] = tgamma(1)*tgamma(k)*(boost::math::digamma(k+1)-boost::math::digamma(k))/tgamma(k+1);
            }
            return mom;
        }

        Vector<T, Dynamic> mmom_jaclog(std::size_t n){
            Vector<T, Dynamic> mmom = Vector<T, Dynamic>::Zero(2*n);
            mmom[0] = tgamma(1)*tgamma(1)*(boost::math::digamma(2)-boost::math::digamma(1))/tgamma(2);
            int sgn = 1;
            for(int i = 1; i < 2*n; i++){
                int k = i+1;
                sgn *= -1;
                mmom[i] = sgn*tgamma(k-1)*tgamma(k)*tgamma(1)/((k)*tgamma(k+1));
            }

            
            return mmom;
        }*/


        Vector<T, Dynamic> mmom_log(std::size_t n){
            Vector<T, Dynamic> mm = Vector<T, Dynamic>::Zero(n);
            double c = 1.0;

            for(int i = 0; i < n; i++){
                if(i == 0){
                    mm[0] = 1;
                }
                else{
                    double p = (i+1)*i;
                    mm[i] = (pow(-1, i)/p)*(pow(tgamma(1), 2));
                }
                mm[i] *=c;
                c *= (0.5*(i+1)/(2*(i+1)-1));
            }

            return mm;
        }


        Matrix<T, Dynamic, Dynamic> chebyshev(std::size_t n, mom, abj){
            // should return an nx2 matrix, with [a,b], with a & b as nx1 vectors
            Matrix<T, Dynamic, Dynamic> abj = Matrix<T, Dynamic, Dynamic>::Zero(n,2);

            ab(0,0) = abj(0,0) + (mom[1]/mom[0]);
            ab(0,1) = mom[0];
            
            // computing sigma
            Matrix<T, Dynamic, Dynamic> sigma = Matrix<T, Dynamic, Dynamic>::Zero(n+1, 2n);

            // initializing first row
            for(int l = 0; l < 2*n; l++){
                sigma(0, l) = mom[l];
            }

            // filling in n following rows
            for(int k = 1; k < n; k++){
                for(int l = k; l < 2*n-k-1; l++){

                    sigma(k,l) = sigma(k-1,l+1) - (ab(k-1,0) - abj(l,0))*sigma(k-1,l)-ab(k-1,1)*sigma(k-2,l) + abj(l,1)*sigma(k-1,l-1);

                }
                // alpha_k
                ab(k, 0) = abj(k, 0)+ (sigma(k,k+1)/sigma(k,k)) - (sigma(k-1,k)/sigma(k-1,k-1));

                // beta_k
                ab(k, 1) = sigma(k,k)/sigma(k-1,k-1);
            }

            return ab;
        }




        /*
        * @method 
        * @brief places the recurrence relation coefficients 'coeffs' in a tridiagonal matrix,
        * following the Golub-Welsch Algorithm
        */
        Matrix<T, Dynamic, Dynamic> tridiagCoeffs(Matrix<T, Dynamic, Dynamic> coeffs, std::size_t n) {
            // argument is a nx2 matrix
            // SIZE CHECK
            assert(coeffs.rows() == n);
            assert(coeffs.cols() == 2);
            
            Matrix<T, Dynamic, Dynamic> tridiag = Matrix<T, Dynamic, Dynamic>::Zero(n, n);

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

        /*
        * @method 
        * @brief computes the nodes & weights of the associated Gauss-Jacobi quadrature rule
        */
        Matrix<T, Dynamic, Dynamic> log_nw(std::size_t n) {

            // finding the coefficients
            double gamma_0 = 1; // unsure, but seems correct, sicne a = b = 0

            //shifted moments
            //Matrix<T, Dynamic, Dynamic> ab = shifted_c_log(2*n);
            //Vector<T, Dynamic> mom = mmom_jaclog(n);
            //Matrix<T, Dynamic, Dynamic> coeffs = chebyshev(n, ab, mom);

            // trying with regular moments
            Vector<T, Dynamic> mom = mom_jaclog(n);
            Matrix<T, Dynamic, Dynamic> coeffs = chebyshev(n, mom);



            // solving the coefficients
            Matrix<T, Dynamic, Dynamic> J_n = tridiagCoeffs(coeffs, n);
            SelfAdjointEigenSolver<Matrix<T, Dynamic, Dynamic>> solve(J_n); // yields much faster computations of high n
            Vector<T, Dynamic> nodes= solve.eigenvalues().real();

            Matrix<T, Dynamic, Dynamic> eigenvecs = solve.eigenvectors().real();

            Vector<T, Dynamic> weights = Vector<T, Dynamic>::Zero(n);

            for(int i = 0; i < n; i++){
                weights[i] = gamma_0*pow(eigenvecs.col(i).normalized()[0], 2);
            }

            Matrix<T, Dynamic, Dynamic> nw = Matrix<T, Dynamic, Dynamic>::Zero(n, 2);
            nw.col(0) = nodes;
            nw.col(1) = weights;
            
            return nw;
            
        }



        template<typename F>
        T operator()(F f) {
            Matrix<T, Dynamic, Dynamic> nw = log_nw(degree);

            T quad = 0;
            for(int i = 0; i < nw.col(0).size(); i++){
                quad+= nw(i, 1)*f(nw(i,0));
            }

            // now the GaussLegendre part

            return quad;
        }

    };







} // Namespace GQLog