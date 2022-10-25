#pragma once

#include "vec.h"
#include "matrix.h"
#include "file_helper.h"
#include <iostream>
#include <map>

extern int num_of_iter;
extern int num_of_calc;
extern int num_of_grad;

double kvad(std::pair<double, double> x);
double roz(std::pair<double, double> x);
vec antigrad_for_kvad(std::pair<double, double> p);
vec antigrad_for_roz(std::pair<double, double> p);
double mzs(double eps, double(*foo)(std::pair<double, double>), std::pair<double, double> point, const vec& antigr, int& schet);
std::pair<std::pair<double, double>, double> method_droblenia_shaga(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_naiscor_spuska(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_sopr_grad_for_kvad(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_sopr_grad_for_roz(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_fletchera_rivsa(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_polaka_rivera(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_newtona(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>), matrix(*matr)(std::pair<double, double> point));
std::pair<std::pair<double, double>, double> method_newtona_naisk(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>), matrix(*matr)(std::pair<double, double> point));
std::pair<std::pair<double, double>, double> method_markvarda(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>), matrix(*matr)(std::pair<double, double> point));
std::pair<std::pair<double, double>, double> method_dfp(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_bfsh(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_paulla(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_pokord_spuska(double(*foo)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_huka_jivsa(double(*foo)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_rozenbroka(double(*foo)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_regular_simplex(double(*foo)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_irregular_simplex(double(*foo)(std::pair<double, double>));
std::pair<std::pair<double, double>, double> method_vnutr_shtraf(double(*foo)(std::pair<double, double>),int kvad_or_roz);
std::pair<std::pair<double, double>, double> method_vnesh_shtraf(double(*foo)(std::pair<double, double>), int kvad_or_roz);


double dot(const vec& antigrad1, const vec& antigrad2);
double dot_with_matrix_Q(const vec& antigrad1, const vec& antigrad2);
double dot_with_matrix_H(const vec& antigrad1, const vec& antigrad2, std::pair<double, double> point);
matrix get_H_for_roz(std::pair<double, double> point);
matrix get_H_for_roz_inverse(std::pair<double, double> point);
matrix get_H_for_kvad(std::pair<double, double> point);
matrix get_H_for_kvad_inverse(std::pair<double, double> point);

double operator-(std::pair<double, double> p1, std::pair<double, double> p2);

vec operator*(matrix A, vec b);


double  vnutr_shtraf_for_kvad1(std::pair<double, double> x);
double  vnutr_shtraf_for_roz1(std::pair<double, double> x);
double  vnutr_shtraf_for_kvad2(std::pair<double, double> x);
double  vnutr_shtraf_for_roz2(std::pair<double, double> x);
double  vnesh_shtraf_for_kvad1(std::pair<double, double> x);
double  vnesh_shtraf_for_roz1(std::pair<double, double> x);
double  vnesh_shtraf_for_kvad2(std::pair<double, double> x);
double  vnesh_shtraf_for_roz2(std::pair<double, double> x);