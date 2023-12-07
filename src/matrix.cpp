#include "matrix.h"
matrix::matrix(std::vector<std::vector<double>> a, int n_)
{
    n = n_;
	A.resize(n);
	for (int i = 0; i < n; ++i)
		A[i].resize(n);
	A = { {a[0][0], a[1][0]},{a[0][1], a[1][1]}};
}
matrix::matrix(base_type_of_matrix type, int size_)
{
	switch (type)
	{
	case base_type_of_matrix::singular:
		matrix::n = size_;
		A.resize(n);
		for (int i = 0; i < n; ++i)
			A[i].resize(n);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				if (i == j) A[i][j] = 1; else A[i][j] = 0;
		break;
	case base_type_of_matrix::zerro:
		matrix::n = size_;
		A.resize(n);
		for (int i = 0; i < n; ++i)
			A[i].resize(n);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				A[i][j] = 0;
		break;
	default:
		matrix::n = size_;
		break;
	}
}
matrix::matrix(const vec& vector)
{
	n = 2;
	A = { {vector.get_x(),0}, {vector.get_y(),0} };
}

std::vector<std::vector<double>> matrix::get_A() const
{
	return A;
}
void matrix::up_tringle()
{

	for (int k = 0; k < A.size() - 1; ++k)
	{
		part_choice(k);
		//	if (A[k][k] == 0) return (false);
		for (int i = k + 1; i < A.size(); ++i)
		{
			double c = A[i][k] / A[k][k];
			for (int j = k; j < A.size() + 1; j++)
			{
				A[i][j] = A[i][j] - c * A[k][j];
			}
		}
	}
	return;
}

void matrix::part_choice(int num_of_current_line)
{
	int num_line_with_min = num_of_current_line;
	for (int i = num_of_current_line + 1; i < A.size(); ++i)
	{
		if (fabs(A[num_of_current_line][num_of_current_line]) < fabs(A[i][num_of_current_line]))
			std::swap(A[num_of_current_line], A[i]);
	}
}

std::vector<double> matrix::get_decision_gaus()
{
	up_tringle();
	matrix B = *this;
	std::vector<double> answ;
	answ.resize(A.size());
	double temp = 0;
	for (int i = A.size() - 1; i >= 0; --i)
	{
		for (int j = i + 1; j < A.size(); ++j)
			temp += A[i][j] * answ[j];
		answ[i] = (A[i][A.size()] - temp) / A[i][i];
		temp = 0;
	}
	*this = B;
	return answ;
}

std::vector<double>  get_decision_gaus_(matrix B, std::vector<double> b)
{
	for (int i = 0; i < B.A.size(); ++i)
		B.A[i].push_back(b[i]);
	return B.get_decision_gaus();
}

matrix matrix::transponate()
{
	matrix B = *this;
	for (int i = 0; i < n; ++i)
		for (int j = i; j < n; ++j)
			std::swap(B.A[i][j], B.A[j][i]);
	return B;
}
matrix matrix::operator*(double a)
{
	matrix B = *this;
	for (int i = 0; i < B.A.size(); ++i)
		for (int j = 0; j < B.A.size(); ++j)
			B.A[i][j] *= a;
	return B;
}
matrix matrix::operator/(double a)
{
	matrix B = *this;
	for (int i = 0; i < B.A.size(); ++i)
		for (int j = 0; j < B.A.size(); ++j)
			B.A[i][j] /= a;
	return B;
}
matrix matrix::operator*(matrix B)
{
	matrix C(base_type_of_matrix::zerro, n);
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			for (int k = 0; k < n; ++k)
				C.A[i][j] += A[i][k] * B.A[k][j];
	return C;
}
matrix matrix::operator+(const matrix& B)
{
	matrix C({ {0,0},{0,0} }, 2);
	for (int i = 0; i < A.size(); ++i)
		for (int j = 0; j < A.size(); ++j)
			C.A[i][j] = A[i][j] + B.A[i][j];
	return C;
}
matrix matrix::operator-(const matrix& B)
{
	matrix C ({{0,0},{0,0}},2);
	for (int i = 0; i < A.size(); ++i)
		for (int j = 0; j < A.size(); ++j)
			C.A[i][j] = A[i][j] - B.A[i][j];
	return C;
}


bool matrix::is_positive()
{
	bool f = ((A[0][0] > 0) && ((A[0][0] * A[1][1] - A[1][0] * A[0][1]) > 0));
	return f;
}