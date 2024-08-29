#include "IPV/range_kutta_method_second_order.hpp"
#include "IPV/euler_method_explict.hpp"
#include "visitor/visitor.hpp"

RangeKuttaSecondOrderMethod::RangeKuttaSecondOrderMethod(double so, double deltat, std::function<double(double)> function, int numbeofstates) : IPV(so, deltat, function, numbeofstates)
{
    this->methodName = "RANGE KUTTA SECOND ORDER";
};

void RangeKuttaSecondOrderMethod::accept(Visitor &visitor) const
{
    visitor.visit(*this);
};

void RangeKuttaSecondOrderMethod::execute()
{
    vector<double> states = {So};
    double Si_bar = 0;
    double Si = 0;
    EulerMethodExplict euler = EulerMethodExplict(So, deltaT, this->functionOfState, 1);
    for (int i = 1; i <= this->numbeOfStates; ++i)
    {
        euler.updateValues(states[i - 1], deltaT, this->functionOfState, 1);
        euler.execute();
        Si_bar = euler.result.back(); // EulerMethods::explicitEulerMethod(states[i - 1], deltaT, function, 1).back();
        Si = states[i - 1] + deltaT * (0.5 * this->functionOfState(states[i - 1]) + 0.5 * this->functionOfState(Si_bar));
        states.push_back(Si);
    }
    this->result = states;
};