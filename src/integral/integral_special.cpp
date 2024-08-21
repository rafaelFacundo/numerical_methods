#include "integral/integral_special.hpp"
#include "visitor/visitor.hpp"

IntegralSpecial::IntegralSpecial(functionWithOneArgument function, double xi, double xf)
    : Xi{xi}, Xf{xf}, functionToIntegrate{function} {};

void IntegralSpecial::accept(Visitor &visitor) const
{
    visitor.visit(*this);
};