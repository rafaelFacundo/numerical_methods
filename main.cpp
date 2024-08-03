#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <limits>
#include <utility>
#include <memory>
#include "include/context_class/context.hpp"
//#include "include/visitor/numerical_method_visitor.hpp"
#include "include/derivate/first_derivate_forward_approach.hpp"
using namespace std;

double functionTeste(double x) {
    return 2.0 * x + 4.0;
} 

int main()
{
    Context teste = Context();
    //NumericalMethodVisitor visi = NumericalMethodVisitor();
    FirstDerivateForwardApproach deri = FirstDerivateForwardApproach(1, 0.1, functionTeste);
    teste.set_strategy(make_unique<FirstDerivateForwardApproach>(1, 0.1, functionTeste));
    teste.callExecute();
   // deri.accept(visi);
    return 0;
}