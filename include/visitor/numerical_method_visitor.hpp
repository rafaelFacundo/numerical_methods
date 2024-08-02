#pragma once
#include <ostream>
#include <iostream>
#include "../derivate/derivate.hpp"
using namespace std;

class NumericalMethodVisitor {
    public:
        /**
         * This method is used to show the result of a derivate 
         *
         * @param derivateToShowResult the derivate class to show the result
         * @return void - this method does not return any value, just prints the result.
         */
        void visit(const Derivate& derivateToShowResult);
};