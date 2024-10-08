#include "integral/integral_special_double_exponential_newton_cotes.hpp"
#include "integral/third_degree_integral_open_newtoncotes.hpp"

#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <math.h>
using namespace std;

IntegralDoubleExponentialNewtonCotes::IntegralDoubleExponentialNewtonCotes(std::function<double(double)> function, double xi, double xf, double tolerance)
    : IntegralSpecial(function, xi, xf), tolerance{tolerance}
{
    this->methodName = "DOUBLE EXPONENTIAL INTEGRAL WITH NEWTON COTES";
};

double IntegralDoubleExponentialNewtonCotes::Xs(double s)
{
    return (this->Xi + this->Xf) / 2 + ((this->Xf - this->Xi) / 2) * tanh((M_PI / 2) * sinh(s));
};

double IntegralDoubleExponentialNewtonCotes::DXs(double s)
{
    return ((this->Xf - this->Xi) / 2) * ((M_PI / 2) * cosh(s)) / pow(cosh((M_PI / 2) * sinh(s)), 2);
};

void IntegralDoubleExponentialNewtonCotes::execute()
{
    auto t = [this](double x)
    {
        return this->functionToIntegrate(this->Xs(x)) * this->DXs(x);
    };

    ThirdDegreeIntegralOpenNewtonCotes integra = ThirdDegreeIntegralOpenNewtonCotes(-3, 6, t, this->tolerance);
    integra.execute();
    this->result = integra.result;
};