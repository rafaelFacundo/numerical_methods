#include "IPV/range_kutta_method_third_order.hpp"
#include "visitor/visitor.hpp"

RangeKuttaThirdOrderMethod::RangeKuttaThirdOrderMethod(double so, double deltat, functionWithOneArgument function, int numbeofstates) : IPV(so, deltat, function, numbeofstates) {};


void RangeKuttaThirdOrderMethod::accept(Visitor& visitor) const {
    visitor.visit(*this);
};

void RangeKuttaThirdOrderMethod::execute() {
   vector<double> states = {So};
    double Si_bar, Si_bar_half, Si = 0;
    for (int i = 1; i <= this->numbeOfStates; ++i)
    {
        Si_bar_half = 0;//EulerMethods::explicitEulerMethod(states[i - 1], deltaT / 2, function, 1).back();
        Si_bar = 0;//EulerMethods::explicitEulerMethod(states[i - 1], deltaT, function, 1).back();
        Si = states[i - 1] + deltaT * ((1.0 / 6.0) * this->functionOfState(states[i - 1]) + (2.0 / 3.0) * this->functionOfState(Si_bar_half) + (1.0 / 6.0) * this->functionOfState(Si_bar));
        states.push_back(Si);
    }
    this->result = states;
};