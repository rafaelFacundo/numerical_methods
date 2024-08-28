#include "derivate/third_derivate_newton_approach.hpp"
#include "visitor/visitor.hpp"
#include <math.h>
using namespace std;

ThirdDerivateNewtonApproach::ThirdDerivateNewtonApproach(double xi, double deltax, std::function<double(double)> function) : Derivate(xi, deltax, function) {};

void ThirdDerivateNewtonApproach::execute()
{
    this->result = (-(5.0 / 2.0) * this->functionToDerive(this->Xi) + 9.0 * this->functionToDerive(this->Xi + this->deltaX) - 12.0 * this->functionToDerive(this->Xi + 2.0 * this->deltaX) + 7.0 * this->functionToDerive(this->Xi + 3.0 * this->deltaX) - (3.0 / 2.0) * this->functionToDerive(this->Xi + 4.0 * this->deltaX)) / pow(this->deltaX, 3);
}

void ThirdDerivateNewtonApproach::accept(Visitor &visitor) const
{
    visitor.visit(*this);
}