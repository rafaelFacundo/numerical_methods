#include "visitor/visitor.hpp"
#include "derivate/derivate.hpp"
#include "integral/integral.hpp"

void Visitor::visit(const Derivate& derivate) {
    cout << "THE DERIVATE WAS CALCULATED FOR: " << '\n';
    cout << "XI      = " << derivate.Xi << '\n';
    cout << "DELTA X = " << derivate.deltaX << '\n';
    cout << "RESLUT  = " << derivate.result << '\n';
}

void Visitor::visit(const Integral& integral) {
  cout << "THE INTEGRAL WAS CALCULATED FOR: " << '\n';
  cout << "XI                    = " << integral.Xi << '\n';
  cout << "DELTA X               = " << integral.deltaX << '\n';
  if(integral.numberOfPartitions != -1) {
    cout << "NUMBER OF PARTITIONS  = " << integral.numberOfPartitions << '\n';
  }else{
    cout << "TOLERANCE             = " << integral.tolerance << '\n';
  }
  cout << "RESULT                = " << integral.result << '\n';
}

void Visitor::teste() {
  cout << "ASKDJLAKSL\n";
}