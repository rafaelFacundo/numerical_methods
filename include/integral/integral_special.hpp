#pragma once
#include "../numerical_method/numerical_method.hpp"
#include "../../types/functions.hpp"

class IntegralSpecial : public NumericalMethod
{
public:
    functionWithOneArgument functionToIntegrate;
    double Xi;
    double Xf;
    double result;

    IntegralSpecial(functionWithOneArgument function, double xi, double xf);

    void accept(Visitor &visitor) const override;
};