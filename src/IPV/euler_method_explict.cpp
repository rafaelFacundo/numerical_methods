#include "IPV/euler_method_explict.hpp"
#include "visitor/visitor.hpp"

EulerMethodExplict::EulerMethodExplict(double so, double deltat, std::function<double(double)> function, int numbeofstates) : IPV(so, deltat, function, numbeofstates)
{
    this->methodName = "EULER EXPLICIT";
};

void EulerMethodExplict::accept(Visitor &visitor) const
{
    visitor.visit(*this);
};

void EulerMethodExplict::execute()
{
    vector<double> states = {this->So};
    double Si = 0;
    for (int i = 1; i <= numbeOfStates; ++i)
    {
        Si = states[i - 1] + deltaT * this->functionOfState(states[i - 1]);
        states.push_back(Si);
    }
    this->result = states;
};