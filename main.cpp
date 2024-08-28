#include <iostream>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <cmath>
#include <limits>
#include <utility>
#include <memory>
#include "include/context_class/context.hpp"
#include "include/derivate/first_derivate_forward_approach.hpp"
#include "include/derivate/first_derivate_backward_approach.hpp"
#include "include/derivate/first_derivate_central_approach.hpp"
#include "include/derivate/second_derivate_forward_approach.hpp"
#include "include/derivate/second_derivate_backward_approach.hpp"
#include "include/derivate/second_derivate_central_approach.hpp"
#include "include/derivate/third_derivate_newton_approach.hpp"

#include "include/visitor/visitor.hpp"

#include "include/integral/integral.hpp"
#include "include/integral/integral_gauss.hpp"
#include "include/integral/integral_special.hpp"
#include "include/integral/first_degree_integral_newtoncotes.hpp"
#include "include/integral/first_degree_integral_open_newtoncotes.hpp"
#include "include/integral/second_degree_integral_open_milnerule.hpp"
#include "include/integral/second_degree_integral_simpson_onethird.hpp"
#include "include/integral/second_degree_integral_gausschebyshev.hpp"
#include "include/integral/second_degree_integral_gausshermite.hpp"
#include "include/integral/second_degree_integral_gausslaguerre.hpp"
#include "include/integral/second_degree_integral_gausslegendre.hpp"
#include "include/integral/third_degree_integral_open_newtoncotes.hpp"
#include "include/integral/third_degree_integral_simpson_threeEighths.hpp"
#include "include/integral/third_degree_integral_gausschebyshev.hpp"
#include "include/integral/third_degree_integral_gausshermite.hpp"
#include "include/integral/third_degree_integral_gausslaguerre.hpp"
#include "include/integral/third_degree_integral_gausslegendre.hpp"
#include "include/integral/integral_special_simple_exponential_newton_cotes.hpp"
#include "include/integral/integral_special_double_exponential_newton_cotes.hpp"

using namespace std;

double sin2x(double x)
{
    return pow(sin(x), 2);
}

double funcao3(double x) {
    return 2 * pow(x, 3) - 3 * pow(x, 2) + 4 * x - 5;
}

double gaussHermite(double x) {
    return  (2 * pow(x,2) + 3 * x + 1);
}

double gaussLaguerre(double x) {
    return  (pow(x,2) + 2 * x + 3);
}

double gaussChebyshev(double x) {
    return (pow(x,2) + 1);
}

double doubleExponential(double x) {
    return 1/sqrt(x);
}

int main()
{

    vector<unique_ptr<NumericalMethod>> testes = {};

  /*   testes.push_back(make_unique<FirstDerivateForwardApproach>(2.0, 0.0001, sin2x));
    testes.push_back(make_unique<FirstDerivateBackwardApproach>(2.0, 0.0001, sin2x));
    testes.push_back(make_unique<FirstDerivateCentralApproach>(2.0, 0.0001, sin2x));

    testes.push_back(make_unique<SecondDerivateCentralApproach>(2.0, 0.0001, sin2x));
    testes.push_back(make_unique<SecondDerivateBackwardApproach>(2.0, 0.0001, sin2x));
    testes.push_back(make_unique<SecondDerivateForwardApproach>(2.0, 0.0001, sin2x));
    testes.push_back(make_unique<ThirdDerivateNewtonApproach>(2.0, 0.0001, sin2x));

     testes.push_back(make_unique<FirstDegreeIntegralNewtonCotes>(1.0, 6.0, sin2x, 10));
     testes.push_back(make_unique<FirstDegreeIntegralOpenNewtonCotes>(1.0, 6.0, sin2x, 10));
     testes.push_back(make_unique<SecondDegreeIntegralOpenMilneRule>(1.0, 6.0, sin2x, 10));

     testes.push_back(make_unique<SecondDegreeIntegralSimpsonOnethird>(1.0, 6.0, sin2x, 10));
     testes.push_back(make_unique<ThirdDegreeIntegralSimpsonThreeEighths>(1.0, 6.0, sin2x, 10));
     testes.push_back(make_unique<ThirdDegreeIntegralOpenNewtonCotes>(1.0, 6.0, sin2x, 10)); 

    testes.push_back(make_unique<SecondDegreeIntegralGuassChebyShev>(1.0, 7.0, gaussChebyshev));
    testes.push_back(make_unique<SecondDegreeIntegralGuassHermite>(1.0, 7.0, gaussHermite));
    testes.push_back(make_unique<SecondDegreeIntegralGuassLaguerre>(1.0, 7.0, gaussLaguerre));
    testes.push_back(make_unique<SecondDegreeIntegralGuassLegendre>(1.0, 7.0, funcao3));  */

    //testes.push_back(make_unique<ThirdDegreeIntegralGuassChebyShev>(1.0, 7.0, gaussChebyshev));
     //testes.push_back(make_unique<ThirdDegreeIntegralGuassHermite>(1.0, 7.0, gaussHermite));
    //testes.push_back(make_unique<ThirdDegreeIntegralGuassLaguerre>(1.0, 7.0, gaussLaguerre));
    /*testes.push_back(make_unique<ThirdDegreeIntegralGuassLegendre>(1.0, 7.0, funcao3));

    /* testes.push_back(make_unique<IntegralSimpleExponentialNewtonCotes>(doubleExponential, 0.0,1.0,0.0001));
    testes.push_back(make_unique<IntegralDoubleExponentialNewtonCotes>(doubleExponential, 0.0,1.0,0.0001)); */







    Context teste = Context();
    Visitor visi = Visitor();
    int i = 0;
    for (auto &method : testes)
    {
        cout << "TESTE " << i << '\n';
        teste.set_strategy(move(method));
        teste.callExecute();
        /* if (teste.nmInstance_) {
            teste.nmInstance_->accept(visi);
        } */
        teste.nmInstance_.get()->accept(visi);
        cout << "=======================\n";
        i += 1;
    }

    /* teste.set_strategy(make_unique<FirstDerivateForwardApproach>());
    teste.callExecute();
    teste.nmInstance_.get()->accept(visi); */
    return 0;
}