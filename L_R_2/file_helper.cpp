#include "file_helper.h"
file_helper::file_helper(void* name_foo1,void* name_foo2)
{
	if (name_foo2 == name_foo1) {
		write_level_lines.open("C:/Users/Артем/Desktop/level_lines_for_kvad.txt");
		write_points.open("C:/Users/Артем/Desktop/points_for_kvad.txt");
		write_simplex.open("C:/Users/Артем/Desktop/simplex_for_kvad.txt");
	}
	else
	{
		write_level_lines.open("C:/Users/Артем/Desktop/level_lines_for_roz.txt");
		write_points.open("C:/Users/Артем/Desktop/points_for_roz.txt");
		write_simplex.open("C:/Users/Артем/Desktop/simplex_for_roz.txt");
	}
}
void file_helper::add_data(std::pair<double, double> x, double y)
{
	write_level_lines << y << std::endl;
	write_points << x.first << " " << x.second << std::endl;
}

void file_helper::add_point(std::pair<double, double> x)
{
	write_points << x.first << " " << x.second << std::endl;
}

void file_helper::add_level_line(double y)
{
	write_level_lines << y << std::endl;
}

void file_helper::add_simplex(const simplex& simp)
{
	auto it = simp.trngl.begin();
	write_simplex << it->second.first << "  " << it->second.second << std::endl;
	it++;
	write_simplex << it->second.first << "  " << it->second.second << std::endl;
	it++;
	write_simplex << it->second.first << "  " << it->second.second << std::endl;
	it = simp.trngl.begin();
	write_simplex << it->second.first << "  " << it->second.second << std::endl;
}

void file_helper::close()
{
	write_level_lines.close();
	write_points.close();
}