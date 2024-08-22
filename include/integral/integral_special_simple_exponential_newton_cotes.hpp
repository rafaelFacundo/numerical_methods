#pragma once
#include "./integral_special.hpp"

class IntegralSimpleExponentialNewtonCotes : public IntegralSpecial {
    public:
        double tolerance;

        IntegralSimpleExponentialNewtonCotes(std::function<double(int)> function, double xi, double xf, double tolerance);

        void execute() override;

        double Xs(double s) override;

        double DXs(double s) override;
};