cmake_minimum_required(VERSION 3.5)

project(opmet)


add_executable(opmet src/file_helper.cpp src/file_helper.h src/matrix.cpp src/matrix.h src/necessary_func.cpp src/necessary_func.h src/simplex.cpp src/simplex.h src/main.cpp src/vec.cpp src/vec.h)

set_target_properties(
      opmet PROPERTIES
      CXX_STANDART 11
      CXX_STANDART_REQUIRED ON)