#pragma once
#include <ostream>
#include <iostream>

using namespace std;

class Derivate;
class Integral;

class Visitor {
    public:
        /**
         * This method is used to show the result of a derivate 
         *
         * @param derivate the derivate class to show thes result
         * @return void - this method does not return any value, just prints the result.
         */
        void visit(const Derivate& derivate);

        /**
         * This method is used to show the result of a integral 
         *
         * @param integral the integral class to show thes result
         * @return void - this method does not return any value, just prints the result.
         */
        void visit(const Integral& integral);

        void teste();
};