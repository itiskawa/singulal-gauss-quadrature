#pragma once
#include <GQJacobi/gauss-rule.hpp>
#include <GQJacobi/gauss-jacobi.hpp>
#include <Eigen/Dense>

using namespace Eigen;

//namespace GQLog{


    /*
    * @author V.B. (alias @itiskawa)
    * @class
    * @brief a class that contains the nodes, weights and number of points of a Gaussian quadrature rule with weight function w(x)=ln(1/x) over ]-1,1[
    * 
    */
    template<typename T>
    class GaussLogRule : public GaussRule<T>{


        protected:

        /*
        * @attributes
        * nodes: associated quadrature rule nodes
        * weights: associated quadrature rule weights
        * degree: number of quadrature points
        */
        std::vector<T> nodes;
        std::vector<T> weights;
        std::size_t degree;

        public:
        /*
        * @method
        * @brief getter for number of quadrature points
        */
        std::size_t getDeg() const{
            return this->degree;
        }

        /*
        * @method
        * @brief getter for quadrature nodes
        */
        std::vector<T> getN() const{
            return this->nodes;
        }

        /*
        * @method
        * @brief getter for quadrature weights
        */
        std::vector<T> getW() const{
            return this->weights;
        }

        /*
        * @constructor
        * @brief default constructor
        */
        GaussLogRule() = default;

        /*
        * @constructor
        * @brief copy constructor
        */
        GaussLogRule(const GaussLogRule& gq){
            this->degree = gq.getDeg();
            for(int i = 0; i < degree; i++){
                this->nodes.push_back(gq.getN()[i]);
                this->weights.push_back(gq.getW()[i]);
            }
        }
        
        /*
        * @constructor 
        * @brief full constructor instantiates the desired quadrature rule by pre-computing the corresponding nodes / weights
        */
        GaussLogRule(std::size_t n) {
            this->degree = n;

            // computing the nodes
            //cout << coeffs.rows() << ", " << coeffs.cols() << endl;
            Matrix<T, Dynamic, Dynamic> nws = nw(n);

            for(int i = 0; i < n; i++){
                nodes.push_back(nws.col(0)[i]);
                weights.push_back(nws.col(1)[i]);
            }
        }
        

        /*
        * @operator
        * @brief template evaluation function
        */
        template<typename F>
        T operator()(F f) {
            T quad = 0;
            // evaluation of integral over ]0,1[, no singularity => use GaussLegendreRule
            GQJacobi::GaussLegendreRule<T> glg(this->degree);
            quad += glg([&](T x){ return log(x+1)*f(x); }, 0, 1);

            // evaluation of integral over ]-1,0[
            for(std::size_t i = 0; i < degree; i++){
                quad -= (weights[i] * std::real(f(nodes[i]-1))) ; // cast to real for cmath functions. Is only meant for f:R->R anyways
            } 
            return quad;
        }

        /*
        * @operator
        * @brief assignment override
        */
        GaussLogRule& operator=(const GaussLogRule& glq){
            if(this != &glq){
                this->degree = glq.getDeg();

                for(size_t i = 0; i < degree; i++){
                    this->nodes.push_back(glq.getN()[i]);
                    this->weights.push_back(glq.getW()[i]);
                }
            }
        }


        protected:
        /*
        * @brief computes the shifted recurrence relation terms
        *
        */
        Matrix<T, Dynamic, Dynamic> shifted_c_log(std::size_t n){
            Matrix<T, Dynamic, Dynamic> abj = Matrix<T, Dynamic, Dynamic>::Zero(n,2);
            Matrix<T, Dynamic, Dynamic> ab = this->c_jacobi(n, 0, 0);
            

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
        Matrix<T, Dynamic, Dynamic> nw(std::size_t n) {

            // finding the coefficients
            double gamma_0 = 1; // given that a=b=0
            Vector<T, Dynamic> mom = mmom_log(2*n);
            Matrix<T, Dynamic, Dynamic> abm = shifted_c_log(2*n);
            Matrix<T, Dynamic, Dynamic> ab = chebyshev(n, mom, abm);
            

            // solving the coefficients
            Matrix<T, Dynamic, Dynamic> J_n = this->tridiagCoeffs(ab, n);
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


    };







//} // Namespace GQLog