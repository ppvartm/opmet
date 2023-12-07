#pragma once
#include <map>
class simplex
{
public:
	std::map<double, std::pair<double, double>> trngl;
	double side_leng;
	double reduct_koef;
	simplex(double(*foo_)(std::pair<double, double>), std::pair<double, double> point);

	void next_step_for_regular(double(*foo_)(std::pair<double, double>), int& k);
	void next_step_for_irregular(double(*foo_)(std::pair<double, double>), int& k);
	double reflect_for_regular(int i, double(*foo)(std::pair<double, double>), int& k);
	double reflect_for_irregular(int i, double(*foo)(std::pair<double, double>), int& k);
	void reduction(double(*foo)(std::pair<double, double>), int& k);
	
	double max_length();

};

