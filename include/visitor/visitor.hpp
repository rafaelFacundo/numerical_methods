#pragma once
#include <ostream>
#include <iostream>

using namespace std;

class Derivate;

class Visitor {
    public:
        /**
         * This method is used to show the result of a derivate 
         *
         * @param derivateToShowResult the derivate class to show thes result
         * @return void - this method does not return any value, just prints the result.
         */
        void visit(const Derivate& derivate);

        void teste();
};