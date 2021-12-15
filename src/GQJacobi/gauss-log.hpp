#pragma once
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cstddef>
#include <type_traits>
using namespace Eigen;

namespace GQLog {

    template<class T>
    class GaussLogRule{


        public:

        std::vector<T> nodes;
        std::vector<T> weights;

        GaussLogRule() = default;

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
            std::cout << "all good" << std::endl;


            abj(0,0) = (1+ab(0,0))/2.;
            abj(0,1) = (ab(0,1))/2;

            for(int i = 1; i < n; i++){
                abj(i,0) = (1+ab(i,0))/2.;
                abj(i,1) = ab(i,1)/4.;
            }
            return abj;
        }


        /*
        * @brief computes modified moments
        *
        */
        Vector<T, Dynamic> mmom(std::size_t n){
            Vector<T, Dynamic> mom = Vector<T, Dynamic>::Zero(n);

            for(int i = 0; i < n; i++){
                mom[i] = pow(-1, i)/((i+1)*i);
            }
            return mom;
        }

        Matrix<T, Dynamic, Dynamic> chebyshev(std::size_t n, Matrix<T, Dynamic, Dynamic> abm, Vector<T, Dynamic> mom){
            // mom must have size 2n
            
            Matrix<T, Dynamic, Dynamic> sig = Matrix<T, Dynamic, Dynamic>::Zero(n+1, 2*n);
            Matrix<T, Dynamic, Dynamic> ab = Matrix<T, Dynamic, Dynamic>::Zero(2*n, 2);

            ab(0,0)=abm(0,0)+mom[1]/mom[0]; 
            ab(0,1)=mom[0];
            std::cout << " so far so good 1" << std::endl;

            sig.col(1) = mom;

            std::cout << " so far so good 2" << std::endl;

            for(int i = 2; i < n+1; i++){
                std::cout << " so far so good i="<<i << std::endl;
                for(int k = i-1; k < 2*n-i+2; k++){
                    sig(i,k)=sig(i-1,k+1)-(ab(i-2,0)-abm(k,0))*sig(n-1,k)-ab(i-2,1)*sig(i-2,k)+abm(k,1)*sig(i-1,k-1);
                }
                ab(i-1,0)=abm(i-1,0)+sig(i,i)/sig(i,i-1)-sig(i-1,i-1)/ sig(i-1,i-2);
                ab(i-1,1)=sig(i,i-1)/sig(i-1,i-2);
            }
            return ab;
        }


        Matrix<T, Dynamic, Dynamic> thang(int n){
            Matrix<T, Dynamic, Dynamic> ab = shifted_c_log(2*n);
            Vector<T, Dynamic> mom = mmom(2*n);
            Matrix<T, Dynamic, Dynamic> nw = chebyshev(n, ab, mom);
            return nw;
        }

        


        
        
        


    };







} // Namespace GQLog