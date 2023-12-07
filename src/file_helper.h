#pragma once
#include <fstream>
#include <utility>
#include "simplex.h"

class file_helper
{
private:
	std::ofstream write_level_lines;
	std::ofstream write_points;
	std::ofstream write_simplex;
public:

	file_helper(void* name_foo1, void* name_foo2);
	void add_data(std::pair<double, double> x, double y);
	void add_point(std::pair<double, double> x);
	void add_level_line( double y);
	void add_simplex(const simplex & simp);
	void close();

};

