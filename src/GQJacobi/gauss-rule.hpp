#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <cstddef>
#include <type_traits>

using namespace Eigen;


template<class T>
class GaussRule {
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
            Vector<T,Dynamic> deg = Vector<T, Dynamic>::LinSpaced(n-1,1, n-1);
            Vector<T,Dynamic> ndeg = (2*deg)+(Vector<T, Dynamic>::Ones(n-1)*(a+b));

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

        double gamma_zero(double a, double b) {
            return (pow(2, a+b+1)*tgamma(a+1)*tgamma(b+1))/(tgamma(a+b+2));
        }

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
        
        // solver for nodes and weights
        Matrix<T, Dynamic, Dynamic> nw(std::size_t n, Matrix<T, Dynamic, Dynamic> ab);
};
