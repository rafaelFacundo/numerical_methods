#include "IPV/adams_bashforth_method_third_order.hpp"
#include "visitor/visitor.hpp"

AdamsBashForthThirdOrder::AdamsBashForthThirdOrder(double so, double deltat, functionWithOneArgument function, int numbeofstates) : IPV(so, deltat, function, numbeofstates) {};


void AdamsBashForthThirdOrder::accept(Visitor& visitor) const {
    visitor.visit(*this);
};

void AdamsBashForthThirdOrder::execute() {
    vector<double> states = {0};//RangeKuttaMethods::rangeKuttaThirdOrderMethod(So, deltaT, function, 2);
    double Si, Si_bar = 0;
    for (int i = 3; i <= this->numbeOfStates; ++i)
    {
        Si_bar = states[i - 1] + (deltaT / 12) * (5 * this->functionOfState(states[i - 3]) - 16 * this->functionOfState(states[i - 2]) + 23 * this->functionOfState(states[i - 1]));
        Si = states[i - 1] + (deltaT / 12) * (8 * this->functionOfState(states[i - 1]) - this->functionOfState(states[i - 2]) + 5 * this->functionOfState(Si_bar));
        states.push_back(Si);
    }
    this->result = states;
};