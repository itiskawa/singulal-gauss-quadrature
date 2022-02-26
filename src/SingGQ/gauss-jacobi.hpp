#pragma once
#include <SingGQ/gauss-rule.hpp>


namespace GQJacobi {

    /**
     * @brief A class that contains the nodes, weights and number of points of a Gauss-Jacobi quadrature rule
     * It inherits from the GaussRule<T> and is also a class template qith parameter T
     * @author V.B. (alias @itiskawa)
     * 
     * @param T : the data type of the quadrature (float, double etc...)
     * @date February 2022
    */
    template <class T>
    class GaussJacobiRule : public GaussRule<T>{

        protected:
        /**
         * @param nodes : nodes of the quadrature rule over ]-1, 1[
         * @param weights : weights of the  quadrature rule over ]-1,1[
         * @param alpha : a parameter of w(x)=(1-x)^a(1+x)^b, must be > -1
         * @param beta : b parameter of w(x)=(1-x)^a(1+x)^b, must be > -1
         * @param degree : number of quadrature points, i.e. size of nodes & weights vector
        */
        std::vector<T> nodes;
        std::vector<T> weights;
        double alpha;
        double beta;
        std::size_t degree;

        public:
        /**
         * @brief getter for number of quadrature points
         * 
         * @return this->degree
        */
        std::size_t getDeg() const{
            return this->degree;
        }

        /**
         * @brief getter for quadrature nodes
         * 
         * @return this->nodes
        */
        std::vector<T> getN() const{
            return this->nodes;
        }

        /**
         * @brief getter for quadrature weights
         * 
         * @return this->weights
        */
        std::vector<T> getW() const{
            return this->weights;
        }

        /**
         * @brief getter for a parameter of w(x)=(1-x)^a(1+x)^b
         * 
         * @return this->alpha
        */
        double getAlpha() const{ return this->alpha; }

        /**
         * @brief getter for b parameter of w(x)=(1-x)^a(1+x)^b
         * 
         * @return this->beta
        */
        double getBeta() const{ return this->beta; }

        /**
         * @brief Default constructor
       */
        GaussJacobiRule() = default;

        /**
         * @brief Copy constructor. 
         * Assigns all fields of gq to this
         * 
         * @param &gq: address of GaussJacobiRule instance to be copied
        */
        GaussJacobiRule(const GaussJacobiRule& gq){
            this->degree = gq.getDeg();
            for(int i = 0; i < degree; i++){
                this->nodes.push_back(gq.getN()[i]);
                this->weights.push_back(gq.getW()[i]);
            }
        }

        /**
         * @brief Constructor
         * Asserts validity of arguments and initializes all fields of this
         * Computes and stores the quadrature nodes and weights
         * 
         * @param n : the desired number of nodes / weightss
         * @param a : the a parameter of the Jacobi weight function
         * @param b : the b parameter of the Jacobi weight function
        */
        GaussJacobiRule(std::size_t n, double a, double b) {

            // validity check
            assert(n>1);
            assert(a > -1);
            assert(b > -1);

            // initializig all fields
            this->degree = n;
            this->alpha = a;
            this->beta = b;
            
            // nodes & weights computation
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nws = nw();

            for(int i = 0; i < n; i++){
                nodes.push_back(nws.col(0)[i]);
                weights.push_back(nws.col(1)[i]);
            }
        }

        /**
         * @brief Applies quadrature rule to a function
         * 
         * @param F : function template type
         * @param f : function template
         * @return numerical value equal to the approximation of I[f,w]
        */
        template<typename F>
        T operator ()(F f) const {
            T quad = 0;
            for(std::size_t i = 0; i < degree; i++){
                // cast to real for cmath functions. Is only meant for f:R->R anyways
                quad += (weights[i] * std::real(f(nodes[i]))) ; 
            } 
            return quad;
        } 

        /**
         * @brief Applies quadrature rule to a function over integration interval ]a,b[
         * Computation is done using affine pullback
         * 
         * @param F : function template type
         * @param f : function template
         * @param a : beginning of integration interval
         * @param b : end of integration interval
         * @return numerical value equal to the approximation of I[f,w] over ]a,b[
        */
        template<typename F>
        T operator ()(F f, T a, T b) const {
            assert(a < b);
            T quad = 0;

            for(std::size_t i = 0; i < degree; i++){
                // affine pullback nodes
                T x_i = 0.5*((1-nodes[i])*a + (1+nodes[i])*b);
                // cast to real for cmath functions. Is only meant for f:R->R anyways
                quad += (weights[i] * std::real(f(x_i))) ; 
            } 
            // affine pullback weight ratio
            quad *= (0.5*(b-a)); 
            return quad;
        }   

        /**
         * @brief Assignment operator 
         * 
         * @param &gq : GaussJacobiRule instance we want to assign from
         * @return *this
        */
        GaussJacobiRule& operator=(const GaussJacobiRule& gq){
            // testing equality first
            if(this != &gq){

                // assigning all fields
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
        /**
         * @brief First moment of the Jacobi weight function
         * 
         * @return the first moment mu_0 of the associated Jacobi weight function
        */
        double gamma_zero() {
            double a = this->alpha;
            double b = this->beta;
            return (pow(2, a+b+1)*tgamma(a+1)*tgamma(b+1))/(tgamma(a+b+2));
            
        }

        /**
          * @brief Computes the nodes and weights of the Gauss-Jacobi quadrature rule
          * The nodes and weights are computed via the Golub-Welsch algorithm using the recurrence relation coefficients
          * of the monic Jaocbi polynomials. The eigenvalue problem is solved using a self-adjoint matrix solver, which
          * significantly speeds up the computation.
          * 
          * @return a (nx2) array with then nodes & weights in the columns
         */
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nw() {

            double gamma_0 = gamma_zero();
            std::size_t n = this->degree;

            /* Computing the recurrence relation coefficients for the tridiagolanisation process
             * Both functions below are from the abstract super-class template GaussRule<T>
            */
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> coeffs = this->c_jacobi(n, this->alpha, this->beta);
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> J_n = this->tridiagCoeffs(coeffs, n);

            // The self-adjoint solver. Significant speedup for large n
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic>> solve(J_n); 
            
            /* Extracting eigenvalues and vectors. Only the real part is considered as 
             * they are guaranteed to be real given the nature of the problem */
            Eigen::Vector<T, Eigen::Dynamic> nodes = solve.eigenvalues().real();
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenvecs = solve.eigenvectors().real();

            // Solving the weights using Golub-Welsch algorithm formula 
            Eigen::Vector<T, Eigen::Dynamic> weights = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
            for(int i = 0; i < n; i++){
                weights[i] = gamma_0*pow(eigenvecs.col(i).normalized()[0], 2);
            }

            // Preparing more compact return type
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nw = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, 2);
            nw.col(0) = nodes;
            nw.col(1) = weights;
            
            return nw;
            
        }

    }; // GaussJacobiRule<T>



    /**
     * @brief A subclass of GaussJacobiRule<T>. Represents an unweighted Gauss quadrature rule i.e. w(x) = 1
     *
    */
    template<typename T>
    class GaussLegendreRule : public GaussJacobiRule<T>{


        public:
        /**
         * @brief Construct a new Gauss Legendre Rule object
         * 
         * @param n : the desired number of quadrature nodes / weights
         */
        GaussLegendreRule(std::size_t n)
        : GaussJacobiRule<T>(n, 0, 0)
        {}

    }; // GaussLegendreRule<T>


    /**
     * @brief A practical subclass of GaussJacobiRule<T>. Represents an Gauss quadrature rule with a Chebyshev weight function.
     * Note that both types of the Chebyshev weight function are accepted, despite the second kind not being singular
     *
    */
    template<typename T>
    class GaussChebyshevRule : public GaussJacobiRule<T>{

        /**
         * @param sgn : The 'sign' of the weight function.
         * -1 will yield the weight function aossiciated to Chebyshev Polynomials of the first kind T(n)
         * whereas +1 will yield the weight function associated to the second kind U(n)
        */
        private:
        int sgn; 

        public:
        /**
         * @brief Constructor of GaussChebyshevRule 
         * 
         * @param n : the desired number of quadrature nodes / weights
         * @param sgn : the sign determining the kind (first or second) of Chebyshev weight function
         */
        GaussChebyshevRule(std::size_t n, int sgn) 
        : GaussJacobiRule<T>(n, sgn*0.5, sgn*0.5)
        {
            // checks that the sign is none other than ±1
            assert(fabs(sgn)==1);
            this->sgn = sgn;
        }

        /**
         * @brief Getter for sgn
         * 
         * @return the sgn field
        */
        int getSgn() const{
            return this->sgn;
        }
    }; // GaussChebyshevRule<T>

} // namespace GQJacobi