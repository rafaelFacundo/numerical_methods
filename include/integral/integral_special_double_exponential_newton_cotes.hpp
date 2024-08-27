#pragma once
#include "./integral_special.hpp"

class IntegralDoubleExponentialNewtonCotes : public IntegralSpecial
{
public:
    double tolerance;

    IntegralDoubleExponentialNewtonCotes(std::function<double(double)> function, double xi, double xf, double tolerance);

    void execute() override;

    double Xs(double s) override;

    double DXs(double s) override;
};