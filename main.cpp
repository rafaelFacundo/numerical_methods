#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <limits>
#include <utility>
#include <memory>
#include "include/context_class/context.hpp"
#include "include/derivate/first_derivate_forward_approach.hpp"
#include "include/visitor/visitor.hpp"
using namespace std;

double functionTeste(double x) {
    return 2.0 * x + 4.0;
} 

int main()
{
    Context teste = Context();
    Visitor visi = Visitor();
    visi.teste();
    FirstDerivateForwardApproach deri = FirstDerivateForwardApproach(1, 0.1, functionTeste);
    teste.set_strategy(make_unique<FirstDerivateForwardApproach>(1, 0.1, functionTeste));
    teste.callExecute();
    teste.nmInstance_.get()->accept(visi);
    return 0;
}