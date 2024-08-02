#include "../../include/visitor/numerical_method_visitor.hpp";

void NumericalMethodVisitor::visit(Derivate& derivateToShowResult) {
    cout << "THE DERIVATE WAS CALCULATED FOR: " << '\n';
    cout << "XI      = " << derivateToShowResult.Xi << '\n';
    cout << "DELTA X = " << derivateToShowResult.deltaX << '\n';
    cout << "RESLUT  = " << derivateToShowResult.result << '\n';
}