#include "visitor/visitor.hpp"
#include "derivate/derivate.hpp"
#include "integral/integral.hpp"
#include "integral/integral_gauss.hpp"
#include "integral/integral_special.hpp"
#include "eigenvalues_eingenvectors/power_methods/power_methods.hpp"
#include "eigenvalues_eingenvectors/house_holder/house_holder.hpp"
#include "eigenvalues_eingenvectors/jacobi/jacobi_method.hpp"
#include "IPV/ipv.hpp"

void Visitor::visit(const Derivate &derivate)
{
  cout << "METHOD " << derivate.methodName << '\n';
  cout << "WAS CALCULATED FOR \n";
  cout << "XI      = " << derivate.Xi << '\n';
  cout << "DELTA X = " << derivate.deltaX << '\n';
  cout << "RESLUT  = " << derivate.result << '\n';
}

void Visitor::visit(const Integral &integral)
{
  cout << "METHOD " << integral.methodName << '\n';
  cout << "WAS CALCULATED FOR \n";
  cout << "XI                    = " << integral.Xi << '\n';
  cout << "DELTA X               = " << integral.deltaX << '\n';
  if (integral.numberOfPartitions != -1)
  {
    cout << "NUMBER OF PARTITIONS  = " << integral.numberOfPartitions << '\n';
  }
  else
  {
    cout << "TOLERANCE             = " << integral.tolerance << '\n';
  }
  cout << "RESULT                = " << integral.result << '\n';
}

void Visitor::visit(const PowerMethods &powerMethod)
{
  /*  cout << "THE POWER METHOD WAS CALCULATED FOR: \n";
   cout << "MATRIX A =\n";
   powerMethod.A.printMatrix();
   cout << "VETOR V0 =\n";
   powerMethod.VectorVo.printMatrix();
   cout << "RESULT   =\n";
   cout << powerMethod.result.eigenValue << '\n'; */
}

void Visitor::visit(const HouseHolderMethod &houseHolder) {
  /* cout << "THE HOUSE HOLDER METHOD WAS CALCULATE FOR: \n";
  houseHolder.A.printMatrix();
  cout << "RESULT       = \n";
  cout << "A HAT MATRIX = \n";
  houseHolder.result.A_hat_matrix.printMatrix();
  cout << "H MATRIX     = \n";
  houseHolder.result.H_matrix.printMatrix(); */
};

void Visitor::visit(const JacobiMethod &jacobi) {
  /* cout << "THE JACOBI METHOD WAS CALCULTATE FOR: \n";
  jacobi.A.printMatrix();
  cout << "REULST eigenvalues = \n";
  jacobi.result.eigenValues.printMatrix(); */
};

void Visitor::visit(const QRMethod &qr) {
  // cout << "THE QR METHOD WAS CALCULATED FOR: \n";
};

void Visitor::visit(const IPV &ipv)
{
  cout << "METHOD " << ipv.methodName << '\n';
  cout << "WAS CALCULATED FOR \n";
  cout << "So      = " << ipv.So << '\n';
  cout << "DELTA T = " << ipv.deltaT << '\n';
  cout << "NUM OF STATES = " << ipv.numbeOfStates << '\n';
  cout << "RESULT VECTOR = [\n";
  for (int i = 0; i < ipv.result.size(); ++i)
  {
    cout << "S{" << i << "} = " << ipv.result[i] << '\n';
  }
  cout << "]\n";
};

void Visitor::visit(const IntegralGauss &integral)
{
  cout << "METHOD " << integral.methodName << '\n';
  cout << "WAS CALCULATED FOR \n";
  cout << "XI                      =  " << integral.Xi << '\n';
  cout << "XF                 =  " << integral.Xf << '\n';
  cout << "ROOTS OF THE INTEGRAL   = [\n";
  for (double root : integral.roots)
  {
    cout << " " << root << "\n";
  }
  cout << "]\n";
  cout << "WEIGHTS OF THE INTEGRAL = [\n";
  for (double weight : integral.weights)
  {
    cout << " " << weight << "\n";
  }
  cout << "]\n";
  cout << "RESULT                  = " << integral.result << '\n';
};

void Visitor::visit(const IntegralSpecial &integral)
{
  cout << "METHOD " << integral.methodName << '\n';
  cout << "WAS CALCULATED FOR \n";
  cout << "XI     =  " << integral.Xi << '\n';
  cout << "XF     =  " << integral.Xf << '\n';
  cout << "RESULT =  " << integral.result << '\n';
};

void Visitor::teste()
{
  cout << "ASKDJLAKSL\n";
}