#include "IPV/adams_bashforth_method_second_order.hpp"
#include "visitor/visitor.hpp"

AdamsBashForthSecondOrder::AdamsBashForthSecondOrder(double so, double deltat, functionWithOneArgument function, int numbeofstates) : IPV(so, deltat, function, numbeofstates) {};


void AdamsBashForthSecondOrder::accept(Visitor& visitor) const {
    visitor.visit(*this);
};

void AdamsBashForthSecondOrder::execute() {
    vector<double> states = {So};
    double Si, Si_bar = 0;
    double S_one = 0;//RangeKuttaMethods::rangeKuttaSecondOrderMethod(So, deltaT, function, 1).back();
    states.push_back(S_one);

    for (int i = 2; i <= this->numbeOfStates; ++i)
    {
        Si_bar = states[i - 1] + (deltaT / 2) * (3 * this->functionOfState(states[i - 1]) - this->functionOfState(states[i - 2]));
        Si = states[i - 1] + deltaT * (0.5 * this->functionOfState(states[i - 1]) + 0.5 * this->functionOfState(Si_bar));
        states.push_back(Si);
    }
    this->result = states;
};