cmake_minimum_required(VERSION 3.10)
project(NumericalMethodsVisitor VERSION 1.0)

include_directories(include)

add_library(${PROJECT_NAME} ./numerical_method_visitor.cpp)