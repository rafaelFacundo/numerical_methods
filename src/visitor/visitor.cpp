#include "visitor/visitor.hpp"
#include "derivate/derivate.hpp"

void Visitor::visit(const Derivate& derivate) {
    cout << "THE DERIVATE WAS CALCULATED FOR: " << '\n';
    cout << "XI      = " << derivate.Xi << '\n';
    cout << "DELTA X = " << derivate.deltaX << '\n';
    cout << "RESLUT  = " << derivate.result << '\n';
}

void Visitor::teste() {
  cout << "ASKDJLAKSL\n";
}