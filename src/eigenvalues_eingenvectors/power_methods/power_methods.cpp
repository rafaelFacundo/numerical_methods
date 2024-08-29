#include "eigenvalues_eingenvectors/power_methods/power_methods.hpp"

PowerMethods::PowerMethods(Matrix a, Matrix vectorvo, double tolerance) : A{a}, VectorVo{vectorvo}, tolerance{tolerance} {}

void PowerMethods::printResult()
{
    cout << "METHOD " << this->methodName << '\n';
    cout << "RESULT = \n";
    cout << "EIGENVALUE = \n";
    cout << this->result.eigenValue << '\n';
    cout << "EIGENVECTOR = \n";
    this->result.eigenVector.printMatrix();
};
