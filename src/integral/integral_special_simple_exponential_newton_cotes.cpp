#include "integral/integral_special_simple_exponential_newton_cotes.hpp"
#include "integral/third_degree_integral_open_newtoncotes.hpp"
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <math.h>
using namespace std;

IntegralSimpleExponentialNewtonCotes::IntegralSimpleExponentialNewtonCotes(std::function<double(double)> function, double xi, double xf, double tolerance)
    : IntegralSpecial(function, xi, xf), tolerance{tolerance}
{
    this->methodName = "SIMPLE EXPONENTIAL INTEGRAL WITH NEWTON COTES";
};

double IntegralSimpleExponentialNewtonCotes::Xs(double s)
{
    return (this->Xi + this->Xf) / 2 + ((this->Xf - this->Xi) / 2) * tanh(s);
};

double IntegralSimpleExponentialNewtonCotes::DXs(double s)
{
    return ((this->Xf - this->Xi) / 2) * (1 / pow(cosh(s), 2));
};

void IntegralSimpleExponentialNewtonCotes::execute()
{
    auto t = [this](double x)
    {
        return this->functionToIntegrate(this->Xs(x)) * this->DXs(x);
    };
    ThirdDegreeIntegralOpenNewtonCotes integra = ThirdDegreeIntegralOpenNewtonCotes(-5, 10, t, this->tolerance);
    integra.execute();
    this->result = integra.result;
};