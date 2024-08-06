#include "IPV/ipv.hpp"

IPV::IPV(double so, double deltat, functionWithOneArgument function, int numbeofstates) : So{so}, deltaT{deltat}, functionOfState{function}, numbeOfStates{numbeofstates} {};