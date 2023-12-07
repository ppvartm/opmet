#pragma once
#include<utility>
#include<vector>

class vec
{
private:
	double x;
	double y;
public:
	vec(double x_, double y_) : x(x_), y(y_) {};
	vec(std::vector<double> vec):x(vec[0]), y(vec[1]) {};

	operator std::vector<double>();

	double norm();
	double get_x() const;
	double get_y() const;

	void set_x(double x_);
	void set_y(double y_);

	vec operator*(double alpha) const;
	std::pair<double, double> operator+(std::pair<double, double> point) const;
	std::pair<double, double> operator-(std::pair<double, double> point) const;
	vec operator+(const vec& vector) const;
	vec operator-(const vec& vector) const;
	friend vec operator*(double alpha, const vec& vector);
	friend std::pair<double, double> operator+(std::pair<double, double> point, const vec& vector);
	friend std::pair<double, double> operator-(std::pair<double, double> point, const vec& vector);
};