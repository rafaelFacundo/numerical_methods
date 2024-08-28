#include "derivate/second_derivate_backward_approach.hpp"
#include "visitor/visitor.hpp"
#include <math.h>
using namespace std;

SecondDerivateBackwardApproach::SecondDerivateBackwardApproach(double xi, double deltax, std::function<double(double)> function) : Derivate(xi, deltax, function) {};

void SecondDerivateBackwardApproach::execute()
{
    this->result = (this->functionToDerive(this->Xi) - 2 * this->functionToDerive(this->Xi - this->deltaX) + this->functionToDerive(this->Xi - 2 * this->deltaX)) / pow(this->deltaX, 2);
}

void SecondDerivateBackwardApproach::accept(Visitor &visitor) const
{
    visitor.visit(*this);
}