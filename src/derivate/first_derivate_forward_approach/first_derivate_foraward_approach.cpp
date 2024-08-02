#include "../../../include/derivate/first_derivate_foward_approach/first_derivate_forward_approach.hpp"

FirstDerivateForwardApproach::FirstDerivateForwardApproach(double xi, double deltax, functionWithOneArgument function) : Derivate(xi, deltax, function) {}


void FirstDerivateForwardApproach::execute() 
{
    this->result = (this->functionToDerive(this->Xi + this->deltaX) - this->functionToDerive(this->Xi)) / this->deltaX;
}

void FirstDerivateForwardApproach::accept(NumericalMethodVisitor visitor) const
{
    visitor.visit(*this);
}