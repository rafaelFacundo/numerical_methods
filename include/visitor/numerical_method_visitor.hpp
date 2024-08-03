#pragma once
#include <ostream>
#include <iostream>
#include "../derivate/derivate.hpp"
#include "../numerical_method/numerical_method.hpp"
using namespace std;

class NumericalMethodVisitor {
    public:
        /**
         * This method is used to show the result of a derivate 
         *
         * @param derivateToShowResult the derivate class to show thes result
         * @return void - this method does not return any value, just prints the result.
         */
        void visit(Derivate& derivateToShowResult);
};