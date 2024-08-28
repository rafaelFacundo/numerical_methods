#include "derivate/first_derivate_backward_approach.hpp"
#include "visitor/visitor.hpp"

FirstDerivateBackwardApproach::FirstDerivateBackwardApproach(double xi, double deltax, std::function<double(double)> function) : Derivate(xi, deltax, function) {};

void FirstDerivateBackwardApproach::execute()
{
    this->result = (this->functionToDerive(this->Xi) - this->functionToDerive(this->Xi - this->deltaX)) / this->deltaX;
}

void FirstDerivateBackwardApproach::accept(Visitor &visitor) const
{
    visitor.visit(*this);
}