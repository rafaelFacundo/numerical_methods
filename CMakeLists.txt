cmake_minimum_required(VERSION 3.10)

# Nome do projeto
project(MyProject)

# Versão do C++ (C++11, C++14, C++17, etc.)
set(CMAKE_CXX_STANDARD 17)

# Diretório de inclusão
include_directories(${CMAKE_SOURCE_DIR}/include)


# Adiciona o executável
add_executable(MyProject main.cpp)
