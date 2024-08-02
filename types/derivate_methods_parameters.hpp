#include "./methods_parameters.hpp";
#include "./functions.hpp";

class DerivadeMethodsParameters : public MethodsParameters {
    public: 
        double Xi;
        double deltaX;
        functionWithOneArgument function;
};