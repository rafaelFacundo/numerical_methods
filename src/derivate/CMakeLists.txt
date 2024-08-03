cmake_minimum_required(VERSION 3.10)
project(Derivates VERSION  1.0)

include_directories(include)

add_library(${PROJECT_NAME} ./derivate.cpp ./first_derivate_foraward_approach.cpp)