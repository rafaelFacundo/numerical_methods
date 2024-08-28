#pragma once
#include "../numerical_method/numerical_method.hpp"
#include "../../types/functions.hpp"
#include <functional>

class Derivate : public NumericalMethod
{
public:
    double deltaX;
    double Xi;
    double result;
    std::function<double(double)> functionToDerive;

    Derivate(double xi, double deltax, std::function<double(double)> function);

    void setXi(double xi);
};