#include "vec.h"

vec::operator std::vector<double>()
{
    std::vector<double>  res(2);
    res[0] = x;
    res[1] = y;
    return res;
}

double vec::norm()
{
    return sqrt(x * x + y * y);
}

vec operator*(double alpha,const vec& vector)
{
   
    return vector*alpha;
}

vec vec::operator*(double alpha) const
{
    
    return vec(alpha * x, alpha * y);
}

std::pair<double, double> vec::operator+(std::pair<double, double> point) const
{
    return std::pair<double, double>(x + point.first, y + point.second);
}

std::pair<double, double> operator+(std::pair<double, double> point, const vec& vector)
{
    return vector + point;
}
vec vec::operator+(const vec& vector) const
{
    return (vec(x + vector.get_x(), y + vector.get_y()));
}
vec vec::operator-(const vec& vector) const
{
    return (vec(x - vector.get_x(), y - vector.get_y()));
}

std::pair<double, double> vec::operator-(std::pair<double, double> point) const
{
    return std::pair<double, double>(x - point.first, y - point.second);
}

std::pair<double, double> operator-(std::pair<double, double> point, const vec& vector)
{
    return { point.first - vector.get_x(),point.second - vector.get_y() };// vector - point;
}
double vec::get_y()const
{
    return y;
}

double vec::get_x()const
{
    return x;
}
void vec::set_x(double x_)
{
    x = x_;
}
void vec::set_y(double y_)
{
    y = y_;
}