#include "../../include/derivate/derivate.hpp"

Derivate::Derivate(double xi, double deltax, std::function<double(double)> function) : Xi{xi}, deltaX{deltax}, functionToDerive{function} {}

void Derivate::setXi(double xi)
{
    this->Xi = xi;
};
