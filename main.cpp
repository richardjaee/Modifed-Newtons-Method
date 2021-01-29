#include <iostream>
#include <cmath>

//euler declaration
const double EulerConstant = std::exp(1.0);

//function declarations
long double function(long double x);
void modifed_newton(long double p0, long double tol, long double max_iter);
long double dfunction(long double x);
long double ddfunction(long double x);

//modified newtons method implementation
void modifed_newton(long double p0, long double tol, long double max_iter){
    int i = 1;
    static long double p_0 = p0;
    while (i < max_iter){
        long double p_i = p_0-((function(p_0)*dfunction(p_0))/((pow(dfunction(p_0),2)-(function(p_0)*ddfunction(p_0)))));
        
        std::cout << 'P' << i << " = " << p_i << std::endl;
        
        long double absolute = abs(p_i-p_0);
        if (absolute < tol){
            std::cout << "Finished in " << i << " iterations" << std::endl;
            return;
        }
        i++;
        p_0 = p_i;
    }
    std::cout << "Method Failed" << std::endl;
}

//function
long double function(long double  x){
    
  // double efun = pow(EulerConstant, x);
  // double xpow = pow(2, -x);
  // double cosin = 2*cos(x);
  // double constant = -6;
  // double squareRoot = sqrt(x);
  // double double variable = 4*x;
  // double result = efun + xpow + cosin + constant; //e^x + 2^-x + 2cosx - 6
  // std::cout << result << std::endl;

  // double cube = pow(x, 3);
  // double double cosin = cos(x);
  // long double result = 1 - 4*x*cos(x) + 2*square + cos(2*x);
    long double square = pow(x, 2);
    long double result = 1-(4*x*cos(x))+(2*square)+cos(2*x);
    return result;
    
}

//second derivative
long double dfunction(long double x){
    
  // double efun = pow(EulerConstant, x);
  // double xpow = pow(2, -x);
  // double cosin = 2*cos(x);
  // double constant = -6;
  // double cube = pow(x, 3);
  // double cosin = cos(x);
  // double squareRoot = sqrt(x);
  // double variable = 0;
  // double result = efun + xpow + cosin + constant; //e^x + 2^-x + 2cosx - 6
  // std::cout << result << std::endl;
    
  // double sinf = sin(x);
  // double square = pow(x, 2);

  // long double result = -2*sin(2*x) + 4*x*sin(x) - 4*cos(x) + 4*x;
    long double result = (-2*sin(2*x))+(4*x*sin(x))-(4*cos(x))+(4*x);
    return result;
}
long double ddfunction(long double x){
    
  // double efun = pow(EulerConstant, x);
  // double xpow = pow(2, -x);
  // double cosin = 2*cos(x);
  // double constant = -6;
  // double cube = pow(x, 3);
  // double cosin = cos(x);
  // double squareRoot = sqrt(x);
  // double variable = 0;
  // double result = efun + xpow + cosin + constant; //e^x + 2^-x + 2cosx - 6
  // std::cout << result << std::endl;
  // double sinf = sin(x);
  // double square = pow(x, 2);

  // long double result = -2*sin(2*x) + 4*x*sin(x) - 4*cos(x) + 4*x;
    long double result = (-4*cos(2*x))+8*sin(x)+(4*x*cos(x))+4;
    return result;
}
int main(int argc, const char * argv[]) {
    //startin p0
    long double p0 = 1.0;
    
    //tolerance
    long double tol = pow(10.0, -5.0);
    
    //iterations
    long double max_iter = 100.0;
    
    //function call
    modifed_newton(p0, tol, max_iter);

}

