#pragma once
#include <GQJacobi/gauss-rule.hpp>


namespace GQJacobi {

    /*
    * @author V.B. (alias @itiskawa)
    * @class
    * @brief a class that contains the nodes, weights and number of points of a Gauss-Jacobi quadrature rule
    * 
    */
    template <class T>
    class GaussJacobiRule : public GaussRule<T>{

        protected:
        /*
        * @attributes
        * nodes: associated quadrature rule nodes
        * weights: associated quadrature rule weights
        * degree: number of quadrature points
        */
        std::vector<T> nodes;
        double alpha;
        double beta;
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
        * @method
        * @brief getter for a parameter of w(x)=(1-x)^a(1+x)^b
        */
        double getAlpha() const{ return this->alpha; }

        /*
        * @method
        * @brief getter for b parameter of w(x)=(1-x)^a(1+x)^b
        */
        double getBeta() const{ return this->beta; }

        /*
        * @constructor
        * @brief default constructor
        */
        GaussJacobiRule() = default;

        /*
        * @constructor
        * @brief copy constructor
        */
        GaussJacobiRule(const GaussJacobiRule& gq){
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
        GaussJacobiRule(std::size_t n, double a, double b) {
            assert(n>1);
            assert(a > -1);
            assert(b > -1);
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nws = nw(n);
            this->degree = n;
            this->alpha = a;
            this->beta = b;

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
        T operator ()(F f) const {
            T quad = 0;
            for(std::size_t i = 0; i < degree; i++){

                quad += (weights[i] * std::real(f(nodes[i]))) ; // cast to real for cmath functions. Is only meant for f:R->R anyways
            } 
            return quad;
        }   

        /*
        * @operator
        * @brief applies affine pullback of the function template f over interval ]a, b[
        *
        */
        template<typename F>
        T operator ()(F f, T a, T b) const {
            T quad = 0;
            for(std::size_t i = 0; i < degree; i++){
                T x_i = 0.5*((1-nodes[i])*a + (1+nodes[i])*b);// affine pullback nodes
                quad += (weights[i] * std::real(f(x_i))) ; // cast to real for cmath functions. Is only meant for f:R->R anyways
            } 
            quad *= (0.5*(b-a)); // affine pullback weight ratio
            return quad;
        }   

        /*
        * @operator
        * @brief assignment override
        */
        GaussJacobiRule& operator=(const GaussJacobiRule& gq){
            if(this != &gq){
                this->degree = gq.getDeg();
                this->alpha = gq.getAlpha();
                this->beta = gq.getBeta();
                
                for(size_t i = 0; i < degree; i++){
                    this->nodes.push_back(gq.getN()[i]);
                    this->weights.push_back(gq.getW()[i]);
                }
            }
            return *this;
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
        * @brief computes the nodes & weights of the associated Gauss-Jacobi quadrature rule
        */
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nw(std::size_t n) {

            double gamma_0 = gamma_zero(this->alpha, this->beta);

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> coeffs = this->c_jacobi(n, this->alpha, this->beta);
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> J_n = this->tridiagCoeffs(coeffs, n);
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>> solve(J_n); // yields much faster computations of high n
            

            Eigen::Vector<T, Eigen::Dynamic> nodes= solve.eigenvalues().real();

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenvecs = solve.eigenvectors().real();

            Eigen::Vector<T, Eigen::Dynamic> weights = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);

            for(int i = 0; i < n; i++){
                weights[i] = gamma_0*pow(eigenvecs.col(i).normalized()[0], 2);
            }

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nw = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, 2);
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
        /*
        * @attributes
        * sgn: -1 will yield the weight function aossiciated to Chebyshev Polynomials of the first kind T(n)
        *      whereas +1 will yield the weight function associated to the second kind U(n)
        */
        private:
        int sgn; 


        public:
        GaussChebyshevRule(std::size_t n, int sgn) 
        : GaussJacobiRule<T>(n, sgn*0.5, sgn*0.5)
        {
            this->sgn = sgn;
        }

        /*
        * @method
        * @brief getter for sign
        */
        int getSgn() const{
            return this->sgn;
        }


    }; // GaussChebyshevRule
    
    


} // namespace GQJacobi