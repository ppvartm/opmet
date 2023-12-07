#include "simplex.h"

simplex::simplex(double(*foo)(std::pair<double, double>), std::pair<double, double> point)
{
	side_leng = 0.5; reduct_koef = 0.5;
	trngl.insert({ foo(point),point });
	std::pair<double, double> p1 = { point.first + ((sqrt(3.) + 1) / (2 * sqrt(2.)) * side_leng),point.second + ((sqrt(3.) - 1) / (2 * sqrt(2.)) * side_leng) };
	std::pair<double, double> p2 = { point.first + ((sqrt(3.) - 1) / (2 * sqrt(2.)) * side_leng),point.second + ((sqrt(3.) + 1) / (2 * sqrt(2.)) * side_leng) };
	trngl.insert({foo(p1),p1});
	trngl.insert({ foo(p2),p2 });
	
	}

void simplex::next_step_for_regular(double(*foo)(std::pair<double, double>), int& k)
{
	double tt = (--trngl.end())->first;
	double new_f = reflect_for_regular(3, foo,k);
	if (new_f < tt)
	{
		trngl.erase(tt);
		return;
	}
	else
	{
		trngl.erase(new_f);
	    tt = (++trngl.begin())->first;
		new_f = reflect_for_regular(2,foo,k);
		if (new_f < tt)
		{
			trngl.erase(tt);
			return;
		}
		else
		{
			trngl.erase(new_f);
			tt = trngl.begin()->first;
			new_f =reflect_for_regular(1, foo,k);
			if (new_f < tt)
			{
				trngl.erase(tt);
				return;
			}
		}
	}
	trngl.erase(new_f);
	reduction(foo,k);
	next_step_for_regular(foo,k);
	return;
}
void simplex::next_step_for_irregular(double(*foo)(std::pair<double, double>), int& k)
{
	double betta = 2.;
	double gamma = 1./ 2;
	double first_top = trngl.begin()->first;
	double second_top = (++trngl.begin())->first;
	double third_top = (--trngl.end())->first;
	double new_top = reflect_for_irregular(3, foo,k);
	if (new_top < first_top)
	{
		trngl.erase(third_top);
		std::pair<double, double> temp_pair = { ((1.-betta)/2)*(trngl.begin()->second.first+(++trngl.begin())->second.first)+betta* trngl.find(new_top)->second.first,
		((1. - betta) / 2)* (trngl.begin()->second.second + (++trngl.begin())->second.second) + betta * trngl.find(new_top)->second.second };
		if (foo(temp_pair) < (++trngl.begin())->first)
		{
			trngl.erase(trngl.begin()->first);
			trngl.insert({ foo(temp_pair), temp_pair });
		}
		return;
	}
	else if (new_top <= second_top)
	{
		trngl.erase(third_top);
		return;
	}
	else if (new_top > second_top)
	{
		std::pair<double, double> temp_pair = {};
		if (third_top >= new_top)
		{
			 temp_pair = { ((1. - gamma) / 2) * (trngl.begin()->second.first + (++trngl.begin())->second.first) + gamma * trngl.find(new_top)->second.first,
		((1. - gamma) / 2) * (trngl.begin()->second.second + (++trngl.begin())->second.second) + gamma * trngl.find(new_top)->second.second };
		}
		else if (third_top < new_top)
		{
			temp_pair = { ((1. - gamma) / 2) * (trngl.begin()->second.first + (++trngl.begin())->second.first) + gamma * trngl.find(third_top)->second.first,
	((1. - gamma) / 2) * (trngl.begin()->second.second + (++trngl.begin())->second.second) + gamma * trngl.find(third_top)->second.second };
		}
		++k;
		if (foo(temp_pair) < foo(trngl.find(third_top)->second))
		{
			trngl.erase(third_top);
			trngl.erase(new_top);
			trngl.insert({ foo(temp_pair), temp_pair });
			return;
		}
	}
	trngl.erase(new_top);
	reduction(foo,k);
	next_step_for_irregular(foo,k);
	return;
}
double simplex::reflect_for_regular(int i, double(*foo)(std::pair<double, double>), int& k)
{
	std::pair<double, double> max_peak;
	auto it = trngl.begin();
	max_peak = (it)->second;
	for (int j = 1; j < i; ++j)
	 max_peak = (++it)->second;

	auto iter_1 = trngl.begin();
	auto iter_2 = trngl.begin();
	{
		if (i == 1)
		{
			iter_1 = ++trngl.begin();
			iter_2 = ++++trngl.begin();
		}
		else if (i == 2)
		{
			iter_1 = trngl.begin();
			iter_2 = ++++trngl.begin();
		}
		else
		{
			iter_1 = trngl.begin();
			iter_2 = ++trngl.begin();
		}
	}
	std::pair<double, double> new_peak = { (iter_1->second.first + iter_2->second.first - max_peak.first),
		(iter_1->second.second + iter_2->second.second - max_peak.second) };
	trngl.insert({ foo(new_peak), new_peak });
	++k;
	return(foo(new_peak));
}
double simplex::reflect_for_irregular(int i, double(*foo)(std::pair<double, double>), int& k)
{
	double aal = 1;
	std::pair<double, double> max_peak;
	auto it = trngl.begin();
	max_peak = (it)->second;
	for (int j = 1; j < i; ++j)
		max_peak = (++it)->second;

	auto iter_1 = trngl.begin();
	auto iter_2 = trngl.begin();
	{
		if (i == 1)
		{
			iter_1 = ++trngl.begin();
			iter_2 = ++++trngl.begin();
		}
		else if (i == 2)
		{
			iter_1 = trngl.begin();
			iter_2 = ++++trngl.begin();
		}
		else
		{
			iter_1 = trngl.begin();
			iter_2 = ++trngl.begin();
		}
	}
	std::pair<double, double> new_peak = { (iter_1->second.first + iter_2->second.first - aal*max_peak.first),
		(iter_1->second.second + iter_2->second.second - aal*max_peak.second) };
	trngl.insert({ foo(new_peak), new_peak });
	++k;
	return(foo(new_peak));
}
void simplex::reduction(double(*foo)(std::pair<double, double>), int& k)
{
	std::pair<double, double> p = (trngl.begin())->second;
	std::pair<double, double> p1 = (++trngl.begin())->second;
	std::pair<double, double> p2 = (++++trngl.begin())->second;
	trngl.erase(foo(p1)); trngl.erase(foo(p2));
	std::pair<double, double> temp1 = { p.first + reduct_koef * (p1.first - p.first),p.second + reduct_koef * (p1.second - p.second) };
	std::pair<double, double> temp2 = { p.first + reduct_koef * (p2.first - p.first),p.second + reduct_koef * (p2.second - p.second) };
	trngl.insert({ foo(temp1), temp1 });
	trngl.insert({foo(temp2), temp2});
	k = k + 2;
	side_leng = side_leng * reduct_koef;
}
double simplex::max_length()
{
	auto p1 = trngl.begin()->second;
	auto p2 = (++trngl.begin())->second;
	auto p3 = (++++trngl.begin())->second;
    
	std::multimap< double,int> ar_of_leng;
	ar_of_leng.insert({ sqrt((p1.first - p2.first) * (p1.first - p2.first) + (p1.second - p2.second) * (p1.second - p2.second)),1 });
	ar_of_leng.insert({ sqrt((p3.first - p2.first) * (p3.first - p2.first) + (p3.second - p2.second) * (p3.second - p3.second)),1 });
	ar_of_leng.insert({ sqrt((p1.first - p3.first) * (p1.first - p3.first) + (p1.second - p3.second) * (p1.second - p3.second)),1 });
	return ar_of_leng.begin()->first;
}