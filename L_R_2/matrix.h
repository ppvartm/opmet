#pragma once
#include <vector>
#include "vec.h"


enum class base_type_of_matrix { singular, zerro };
enum class type_of_matrix_from_file { square, n_slau };

class matrix
{
private:
	int n;
	std::vector<std::vector<double>> A;

public:
	matrix(std::vector<std::vector<double>> a, int n_);
	matrix(base_type_of_matrix type, int size_);
	matrix(const vec& vector);

	std::vector<double> get_decision_gaus();
	std::vector<double> friend get_decision_gaus_(matrix B, std::vector<double> b);
	std::vector<std::vector<double>> get_A() const;
	void up_tringle();
	void part_choice(int num_of_current_line);
	
	matrix transponate();
	matrix operator*(double a);
	matrix operator/(double a);
	matrix operator*(matrix B);
	matrix operator+(const matrix& B);
	matrix operator-(const matrix& B);
	bool is_positive();
};

