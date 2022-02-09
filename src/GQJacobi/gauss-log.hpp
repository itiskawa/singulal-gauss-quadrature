#pragma once
#include <GQJacobi/gauss-rule.hpp>
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cstddef>
#include <type_traits>
using namespace Eigen;

//namespace GQLog{

    template<typename T>
    class GaussLogRule{


        public:

        std::vector<T> nodes;
        std::vector<T> weights;
        std::size_t degree;

        GaussLogRule() = default;

        GaussLogRule(std::size_t n) {
            this->degree = n;

            // computing the nodes
            Vector<T, Dynamic> mom = mmom_log(2*n);
            Matrix<T, Dynamic, Dynamic> abm = shifted_c_log(2*n);
            Matrix<T, Dynamic, Dynamic> coeffs = chebyshev(n, mom, abm);
            cout << coeffs.rows() << ", " << coeffs.cols() << endl;
            Matrix<T, Dynamic, Dynamic> nw = nw(n, coeffs);

            for(int i = 0; i < n; i++){
                nodes.push_back(nw.col(0)[i]);
                weights.push_back(nw.col(1)[i]);
            }
        }

        /*
        * @method 
        * @brief computes the recurrence relation coefficents (alpha_n, beta_n) of the
        * monic polynomials associated to the Jacobi weight function
        */
       /*
        
        /*
        * @brief computes the shifted recurrence relation terms
        *
        */
        Matrix<T, Dynamic, Dynamic> shifted_c_log(std::size_t n){
            Matrix<T, Dynamic, Dynamic> abj = Matrix<T, Dynamic, Dynamic>::Zero(n,2);
            Matrix<T, Dynamic, Dynamic> ab = c_jacobi<T>(n, 0, 0);


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

        /*
        * @method
        * @brief returns the recurrence coefficients of ln(1/t) for t in ]0,1[ as a nx2 array
        *
        */
        Matrix<T, Dynamic, Dynamic> chebyshev(std::size_t n,Vector<T, Dynamic> mom,Matrix<T, Dynamic, Dynamic> abj){
            // should return an nx2 matrix, with [a,b], with a & b as nx1 vectors
            Matrix<T, Dynamic, Dynamic> ab = Matrix<T, Dynamic, Dynamic>::Zero(n,2);

            ab(0,0) = abj(0,0) + (mom[1]/mom[0]);
            ab(0,1) = mom[0];

            T s_1 = 0.;
            
            // computing sigma
            Matrix<T, Dynamic, Dynamic> sigma = Matrix<T, Dynamic, Dynamic>::Zero(n+1, 2*n);

            // initializing first row
            for(int l = 0; l < 2*n; l++){
                sigma(0, l) = mom[l];
            }

            // filling in n following rows
            for(int k = 1; k < n; k++){
                for(int l = k; l < 2*n-k; l++){
                    if(k > 1){ s_1 = sigma(k-2, l); }
                    sigma(k,l) = sigma(k-1,l+1) - (ab(k-1,0) - abj(l,0))*sigma(k-1,l)-ab(k-1,1)*s_1 + abj(l,1)*sigma(k-1,l-1);

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
        * @brief computes the nodes & weights of the associated Gauss-Jacobi quadrature rule
        */
        Matrix<T, Dynamic, Dynamic> nw(std::size_t n, Matrix<T, Dynamic, Dynamic> ab) {

            // finding the coefficients
            double gamma_0 = 1; // given that a=b=0

            
            

            // solving the coefficients
            Matrix<T, Dynamic, Dynamic> J_n = tridiagCoeffs(ab, n);
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
            T quad = 0;
            for(std::size_t i = 0; i < degree; i++){
                quad += (weights[i] * std::real(f(nodes[i]))) ; // cast to real for cmath functions. Is only meant for f:R->R anyways
            } 
            // now the GaussLegendre part

            return quad;
        }

    };







//} // Namespace GQLog