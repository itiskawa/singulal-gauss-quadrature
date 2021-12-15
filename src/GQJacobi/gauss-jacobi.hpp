#pragma once
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cstddef>
#include <type_traits>
using namespace Eigen;

namespace GQJacobi{

    /*
    * @author V.B (@itiskawa)
    * @class
    * @brief a class that contains the nodes, weights and number of points of a Gauss-Jacobi quadrature rule
    * 
    */
    template <class T>
    class GaussJacobiRule{

        
        public:
        /*
        * @attributes
        * nodes: associated quadrature rule nodes
        * weights: associated quadrature rule weights
        * degree: number of quadrature points
        */
        std::vector<T> nodes;
        std::vector<T> weights;
        std::size_t degree;


        /*
        * @constructor
        * @brief default constructor
        */
        GaussJacobiRule() = default;

        /*
        * @constructor
        * @brief copy constructor
        */
        GaussJacobiRule(const GaussJacobiRule& g){
            this->degree = g.degree;
            for(int i = 0; i < degree; i++){
                this->nodes[i] = g.nodes[i];
                this->weights[i] = g.weights[i];
            }
        }

        /*
        * @constructor 
        * @brief full constructor instantiates the desired quadrature rule by pre-computing the corresponding nodes / weights
        */
        GaussJacobiRule(std::size_t n, double a, double b) {
            assert(n>1);
            assert(a > -1);
            assert(b > -1);
            Matrix<T, Dynamic, Dynamic> nw = jacobi_nw(n, a, b);
            this->degree = n;

            for(int i = 0; i < n; i++){
                nodes.push_back(nw.col(0)[i]);
                weights.push_back(nw.col(1)[i]);
            }
        }



        /*
        * @operator
        * @brief template evaluation function
        */
        template<typename F>
        T operator ()(F f) const {
            T quad = 0;
            for(std::size_t i = 0; i < degree; i++){

                quad += (weights[i] * std::real(f(nodes[i]))) ; // cast to real for cmath functions. Is only meant for f:R->R anyways
            } 
            return quad;
        }   

        /*
        * @operator
        * @brief assignment override
        */
        GaussJacobiRule& operator=(const GaussJacobiRule& gq){
            if(this != &gq){
                this->degree = gq.degree;
                for(size_t i = 0; i < degree; i++){
                    this->nodes.push_back(gq.nodes[i]);
                    this->weights.push_back(gq.weights[i]);
                }
            }

        }



        protected:
        /*
        * @method
        * @brief first moment of the Jacobi weight function
        */
        double gamma_zero(double a, double b) {
            return (pow(2, a+b+1)*tgamma(a+1)*tgamma(b+1))/(tgamma(a+b+2));
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
        Matrix<T, Dynamic, Dynamic> jacobi_nw(std::size_t n, double a, double b) {

            double gamma_0 = gamma_zero(a, b);

            Matrix<T, Dynamic, Dynamic> coeffs = c_jacobi(n, a, b);
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

    
    }; // GaussJacobiRule



    /*
    * @class
    * @brief a subclass of GaussJacobiRule. Represents an unweighted Gauss quadrature rule
    *
    */
    template<typename T>
    class GaussLegendreRule : public GaussJacobiRule<T>{


        public:
        GaussLegendreRule(std::size_t n)
        : GaussJacobiRule<T>(n, 0, 0)
        {}

    }; // GaussLegendreRule


    /*
    * @class
    * @brief a subclass of GaussJacobiRule. Represents an Gauss quadrature rule with a Chebyshev weight function.
    * sgn = -1 or +1. -1 yields polynomials of first kind T(n). +1 yields the second kind U(n)
    *
    */
    template<typename T>
    class GaussChebyshevRule : public GaussJacobiRule<T>{
        public:
        GaussChebyshevRule(std::size_t n, int sgn) 
        : GaussJacobiRule<T>(n, sgn*0.5, sgn*0.5)
        {
            this->sgn = sgn;
        }

        /*
        * @attributes
        * sgn: -1 will yield the weight function aossiciated to Chebyshev Polynomials of the first kind T(n)
        *      whereas +1 will yield the weight function associated to the second kind U(n)
        */
        public:
        int sgn; 

    }; // GaussChebyshevRule
    
    


} // namespace GQJacobi cool