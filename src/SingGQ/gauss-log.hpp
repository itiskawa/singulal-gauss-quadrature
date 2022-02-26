#pragma once
#include <SingGQ/gauss-rule.hpp>
#include <SingGQ/gauss-jacobi.hpp>


namespace GQLog{

    /**
     * @brief A class that contains the nodes, weights and number of points of a Gaussian quadrature rule with weight function w(x)=ln(1/x) over ]-1,1[
     * It inherits from the GaussRule<T> and is also a class template qith parameter T
     * @author V.B. (alias @itiskawa)
     * 
     * @param T : the data type of the quadrature (float, double etc...)
     * @date February 2022
    */
    template<typename T>
    class GaussLogRule : public GaussRule<T>{


        protected:

        /**
         * @param nodes : nodes of the quadrature rule Q[., ln(1/x)] over ]0,1[
         * @param weights : weights of the quadrature Q[., ln(1/x)] rule over ]0,1[
         * @param degree : number of quadrature points, i.e. size of nodes & weights vector
        */
        std::vector<T> nodes;
        std::vector<T> weights;
        std::size_t degree;

        public:
        /**
         * @brief Getter for number of quadrature points
         * 
         * @return this->degree
         */
        std::size_t getDeg() const{
            return this->degree;
        }

        /**
         * @brief Getter for quadrature nodes
         * 
         * @return this->nodes
         */
        std::vector<T> getN() const{
            return this->nodes;
        }

        /**
         * @brief Getter for quadrature weights
         * 
         * @return this->weights
         */
        std::vector<T> getW() const{
            return this->weights;
        }

        /**
         * @brief Default constructor for GaussLogRule 
         * 
         */
        GaussLogRule() = default;

        /**
         * @brief Copy constructor. 
         * Assigns all fields of gql to this
         * 
         * @param &gql: address of GaussLogRule instance to be copied
        */
        GaussLogRule(const GaussLogRule& gql){
            this->degree = gql.getDeg();
            for(int i = 0; i < degree; i++){
                this->nodes.push_back(gql.getN()[i]);
                this->weights.push_back(gql.getW()[i]);
            }
        }
        
        /**
         * @brief Constructor of a new GaussLogRule object
         * 
         * @param n : desired number of quadrature nodes / weights
         */
        GaussLogRule(std::size_t n) {
            this->degree = n;

            // Computing the nodes
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nws = nw(n);

            for(int i = 0; i < n; i++){
                nodes.push_back(nws.col(0)[i]);
                weights.push_back(nws.col(1)[i]);
            }
        }
        

        /**
         * @brief Applies quadrature rule to a function
         * This evaluation is a bit tricker, as the desired weight function is actually ln(x+1) over ]-1,1[, not ln(1/x) over ]0,1[
         * Because of this, half of the quadrature is solved via a GaussLegendreRule, and the other half by this
         * 
         * @param F : function template type
         * @param f : function template
         * @return numerical value equal to the approximation of I[f,ln(x+1)]
        */
        template<typename F>
        T operator()(F f) {
            T quad = 0;
            
            // evaluation of integral over ]0,1[, no singularity => use GaussLegendreRule
            GQJacobi::GaussLegendreRule<T> glg(this->degree);
            quad += glg([&](T x){ return log(x+1)*f(x); }, 0, 1);


            // evaluation of integral over ]-1,0[
            for(std::size_t i = 0; i < degree; i++){
                // cast to real for cmath functions. Is only meant for f:R->R anyways
                quad -= (weights[i] * std::real(f(nodes[i]-1))) ; 
            } 
            return quad;
        }


        template<typename F>
        T operator()(F f, T a, T b) {
            assert(a < b);
            assert(a >= -1);
            T quad = 0;
            
            // evaluation of integral over ]0,1[, no singularity => use GaussLegendreRule
            GQJacobi::GaussLegendreRule<T> glg(this->degree);
            quad += glg([&](T x){ return log(x+1)*f(x); }, 0, 1);


            // evaluation of integral over ]-1,0[
            for(std::size_t i = 0; i < degree; i++){
                // cast to real for cmath functions. Is only meant for f:R->R anyways
                quad -= (weights[i] * std::real(f(nodes[i]-1))) ; 
            } 
            return quad;
        }



        /**
         * @brief Assignment operator 
         * 
         * @param &glq : GaussLogiRule instance we want to assign from
         * @return *this
        */
        GaussLogRule& operator=(const GaussLogRule& glq){
            if(this != &glq){
                this->degree = glq.getDeg();

                for(size_t i = 0; i < this->degree; i++){
                    this->nodes.push_back(glq.getN()[i]);
                    this->weights.push_back(glq.getW()[i]);
                }
            }
            return *this;
        }


        protected:
        /**
         * @brief Computes the shifted recurrence relation coefficients of the Jacobi weight function
         * The computation follows W. Gautschi's r_jacobi01.m method, available at: https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
         * 
         * @param n : desired number of coefficients
         * @return a (nx2) array with the shifted monic Jacobi recurrence relation coefficients
         */
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> shifted_c_log(std::size_t n){
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> abj = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n,2);

            // Computing the regular recurrence relation coefficients
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ab = this->c_jacobi(n, 0, 0);
            
            // Applying W. Gautschi's method
            for(int i = 0; i < n; i++){
                abj(i,0) = (1+ab(i,0))/2.;
            }
            abj(0, 1) = ab(0,1)/2.;
            for(int i = 1; i < n; i++){
                abj(i,1) = ab(i,1)/4.;
            }
            return abj;
        }

        /**
         * @brief Computes the modified moments (mixed moments) of ln(1/x) with respect to the Jacobi Polynomials
         * Routine is found in mmom_log.m by W. Gautschi, available at: https://www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
         * 
         * @param n : the desired number of moments
         * @return The n first modified moments
         */
        Eigen::Vector<T, Eigen::Dynamic> mmom_log(std::size_t n){
            Eigen::Vector<T, Eigen::Dynamic> mm = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
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

        /**
         * @brief The modified Chebyshev Algorithm. It computes the recurrence relation coefficients of the orthogonal polynomial w.r. to w(x) = ln(1/x)
         * The mathematics of the computation are explained in Chapter 2 of W. Gautschi's 'Orthogonal Polynomials: Computation and Approximation'
         * 
         * 
         * @param n : the desired amount of recurrence relation coefficients
         * @param mom : the 2n first modified moments
         * @param abj : the 2n first shifted recurrence relation coefficients
         * @return The recurrence relation coefficients of the orthogonal polynomial w.r. to w(x) = ln(1/x) 
         */
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> chebyshev(std::size_t n,Eigen::Vector<T, Eigen::Dynamic> mom, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> abj){
            // size check
            assert(2*n == mom.rows());
            assert(2*n == abj.rows());

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ab = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n,2);

            ab(0,0) = abj(0,0) + (mom[1]/mom[0]);
            ab(0,1) = mom[0];

            T s_1 = 0.;
            
            // Computing sigma-matrix
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> sigma = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n+1, 2*n);

            // Initializing first row
            for(int l = 0; l < 2*n; l++){
                sigma(0, l) = mom[l];
            }

            // Filling in n following rows
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

        /**
         * @brief computes the nodes & weights of the associated Gauss-Jacobi quadrature rule
         * 
         * @param n : desired number of quadrature nodes / weights
        */
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nw(std::size_t n) {

            // Finding the coefficients
            double gamma_0 = 1; // given that a=b=0
            Eigen::Vector<T, Eigen::Dynamic> mom = mmom_log(2*n);
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> abm = shifted_c_log(2*n);
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> ab = chebyshev(n, mom, abm);
            

            // Solving for the coefficients with Golub-Welsch
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> J_n = this->tridiagCoeffs(ab, n);
           
            // The self-adjoint solver. Significant speedup for large n
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> solve(J_n); // yields much faster computations of high n

            /* Extracting eigenvalues and vectors. Only the real part is considered as 
             * they are guaranteed to be real given the nature of the problem */
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> eigenvecs = solve.eigenvectors().real();
            Eigen::Vector<T, Eigen::Dynamic> nodes= solve.eigenvalues().real();

            // Solving the weights using Golub-Welsch algorithm formula 
            Eigen::Vector<T, Eigen::Dynamic> weights = Eigen::Vector<T, Eigen::Dynamic>::Zero(n);
            for(int i = 0; i < n; i++){
                weights[i] = gamma_0*pow(eigenvecs.col(i).normalized()[0], 2);
            }

            // Preparing return type 
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> nw = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(n, 2);

            nw.col(0) = nodes;
            nw.col(1) = weights;
            
            return nw;
        }
    }; // GaussLogRule<T>


} // Namespace GQLog