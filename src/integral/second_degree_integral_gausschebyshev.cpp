#include "integral/second_degree_integral_gausschebyshev.hpp"
#include "visitor/visitor.hpp"
#include <math.h>
using namespace std;

SecondDegreeIntegralGuassChebyShev::SecondDegreeIntegralGuassChebyShev(double xi, double xf, std::function<double(double)> function) 
: IntegralGauss(xi, xf, function, vector<double>{-(1/sqrt(2)), (1/sqrt(2))}, vector<double>{(M_PIf/2), (M_PIf/2)}) {};

/**
 * This method is the accept method for the visitor pattern
 * used to get the result of the calculation of a method
 *
 *
 * @return void - this method does not return any value, just prints the result.
 */
void SecondDegreeIntegralGuassChebyShev::accept(Visitor &visitor) const
{
    visitor.visit(*this);
};

double SecondDegreeIntegralGuassChebyShev::getXk(int rootIndex)
{
    return this->roots[rootIndex];
};

/**
 * Used by the context class to c\all the execution of some numerical method
 * This method does not receive any parameter
 *
 * @return void - This method does not return any value, just sets the result on the result attribute
 */
void SecondDegreeIntegralGuassChebyShev::execute()
{
    double sum = 0;
    for (int i = 0; i < roots.size(); ++i)
    {
        sum += this->functionToIntegrate(this->getXk(i)) * this->weights[i];
    }
    this->result = sum;
};