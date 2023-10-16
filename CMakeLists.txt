cmake_minimum_required(VERSION 3.5)

project(opmet)


add_executable(opmet L_R_2/file_helper.cpp L_R_2/file_helper.h L_R_2/matrix.cpp L_R_2/matrix.h L_R_2/necessary_func.cpp L_R_2/necessary_func.h L_R_2/simplex.cpp L_R_2/simplex.h L_R_2/main.cpp L_R_2/vec.cpp L_R_2/vec.h)

set_target_properties(
      opmet PROPERTIES
      CXX_STANDART 11
      CXX_STANDART_REQUIRED ON)