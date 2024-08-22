#pragma once
#include "../numerical_method/numerical_method.hpp"

#include <functional>

class IntegralSpecial : public NumericalMethod
{
public:
    std::function<double(int)> functionToIntegrate;
    double Xi;
    double Xf;
    double result;

    IntegralSpecial(std::function<double(int)> function, double xi, double xf);

    void accept(Visitor &visitor) const override;

    virtual double Xs(double s) = 0;

    virtual double DXs(double s) = 0;
};