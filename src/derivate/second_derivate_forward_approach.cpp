#include "derivate/second_derivate_forward_approach.hpp"
#include "visitor/visitor.hpp"
#include <math.h>
using namespace std;

SecondDerivateForwardApproach::SecondDerivateForwardApproach(double xi, double deltax, functionWithOneArgument function) : Derivate(xi,deltax, function) {};

void SecondDerivateForwardApproach::execute() 
{
    this->result = (this->functionToDerive(this->Xi + 2.0 * this->deltaX) - 2.0 * this->functionToDerive(this->Xi + this->deltaX) + this->functionToDerive(this->Xi)) / pow(this->deltaX, 2);
}

void SecondDerivateForwardApproach::accept(Visitor& visitor) const
{
    visitor.visit(*this);
}