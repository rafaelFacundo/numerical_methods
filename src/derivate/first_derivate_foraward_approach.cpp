#include "derivate/first_derivate_forward_approach.hpp"
#include "visitor/visitor.hpp"
#include <iostream>
#include <ostream>
using namespace std;

FirstDerivateForwardApproach::FirstDerivateForwardApproach(double xi, double deltax, std::function<double(double)> function) : Derivate(xi, deltax, function)
{
    this->methodName = "FIRST DERIVATE FORWARD APPROACH";
}

void FirstDerivateForwardApproach::execute()
{
    this->result = (this->functionToDerive(this->Xi + this->deltaX) - this->functionToDerive(this->Xi)) / this->deltaX;
}

void FirstDerivateForwardApproach::accept(Visitor &visitor) const
{
    visitor.visit(*this);
}