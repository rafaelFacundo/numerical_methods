#pragma once
#include "../../types/functions.hpp"
#include "../numerical_method/numerical_method.hpp"
#include <vector>
#include <functional>
using namespace std;

class IPV : public NumericalMethod
{
public:
    double So;
    double deltaT;
    int numbeOfStates;
    std::function<double(double)> functionOfState;
    vector<double> result;

    IPV(double so, double deltat, std::function<double(double)> function, int numbeofstates);

    void updateValues(double so, double deltat, std::function<double(double)> function, int numbeofstates);
};