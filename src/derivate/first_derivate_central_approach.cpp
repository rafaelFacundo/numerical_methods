#include "derivate/first_derivate_central_approach.hpp"
#include "visitor/visitor.hpp"

FirstDerivateCentralApproach::FirstDerivateCentralApproach(double xi, double deltax, std::function<double(double)> function) : Derivate(xi, deltax, function)
{
    this->methodName = "FIRST DERIVATE CENTRAL APPROACH";
};

void FirstDerivateCentralApproach::execute()
{
    this->result = (this->functionToDerive(this->Xi + this->deltaX) - this->functionToDerive(this->Xi - this->deltaX)) / (2 * this->deltaX);
}

void FirstDerivateCentralApproach::accept(Visitor &visitor) const
{
    visitor.visit(*this);
}