#include "necessary_func.h"

extern const int alpha ;
extern const int n;
extern int num_of_iter ;
extern int num_of_calc ;
extern int num_of_grad;
extern int num_of_second_d;
extern double r_for_vnutr;
extern double r_for_vnesh;
extern std::string name_of_current_foo;
extern  const  double kappa0;
extern const double omega;
extern const double lambda0;
extern const double nu;
extern const double eps;
extern  const std::pair<double, double> nachalnaya_tochka ;
extern const matrix Q;


double kvad(std::pair<double, double> p)
{
	return 2 * p.first * p.first - 4 * p.first * p.second + 5 * p.second * p.second - 4 * sqrt(5) * (p.first - p.second) + 4;
	//return 6 * p.first * p.first - 4 * p.first * p.second + 3 * p.second * p.second + 4 * sqrt(5) * (p.first + 2 * p.second) + 22;
}
double roz(std::pair<double, double> p)
{
	return alpha * (p.first * p.first - p.second) * (p.first * p.first - p.second) + (p.first - 1) * (p.first - 1);
}
vec antigrad_for_kvad(std::pair<double, double> p)
{
	return vec(-(4 * p.first - 4 * p.second - 4 * sqrt(5)), -(-4 * p.first + 10 * p.second + 4 * sqrt(5)));
}
vec antigrad_for_roz(std::pair<double, double> p)
{
	return vec(-(2 * alpha * (p.first * p.first - p.second) * 2 * p.first + 2 * (p.first - 1)), -(-2 * alpha * (p.first * p.first - p.second)));
}
double mzs(double eps, double(*foo)(std::pair<double, double>), std::pair<double, double> point, const vec& antigr, int& schet)
{
	double tao = (1 + sqrt(5)) / 2;
	double left = -5;
	double right =5;
	double x1 = left + ((right - left) - (right - left) / tao);
	double x2 = left + (right - left) / tao;
	double foo1 = foo(point + x1 * antigr);
	double foo2 = foo(point + x2 * antigr);
	++schet;
	++schet;
	while (right - left > eps)
	{
		if (foo1 < foo2)
		{
			right = x2;
			x2 = x1;
			x1 = left + right - x2;
			foo2 = foo1;
			foo1 = foo(point + x1 * antigr);
			++schet;
		}
		else
		{
			left = x1;
			x1 = x2;
			x2 = left + right - x1;
			foo1 = foo2;
			foo2 = foo(point + x2 * antigr);
			++schet;
		}


	}
	return (right + left) / 2;
}
std::pair<std::pair<double, double>, double> method_droblenia_shaga(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));//открытие файла
	std::pair<double, double> x0 = nachalnaya_tochka;
	double f1 = foo(x0);
	write.add_data(x0, f1); //записываем начальную точку в файл
	double f2 = 0.;
	std::pair<double, double> x = { };

	double kappa = kappa0;
	while (antigrad(x0).norm() > eps)
	{
		x = x0 + kappa * antigrad(x0);
		f2 = foo(x);
		++num_of_calc;
		while (f1 - f2 < omega * kappa * antigrad(x0).norm() * antigrad(x0).norm())
		{
			kappa = nu * kappa;
			x = x0 + kappa * antigrad(x0);
			f2 = foo(x);
			++num_of_calc;
		}
		kappa = kappa0;
		x0 = x;
		f1 = f2;
		write.add_data(x0, f1);
		++num_of_iter;
		++num_of_grad;
	}
	write.close();
	return { {x0.first,x0.second}, f1 };
}
std::pair<std::pair<double, double>, double> method_naiscor_spuska(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	double f1 = foo(x0);
	write.add_data(x0, f1);
	std::pair<double, double> x = { };
	double kappa = kappa0;
	while (antigrad(x0).norm() > eps)
	{
		kappa = mzs(eps / 1000, foo, x0, antigrad(x0), num_of_calc);
		//foo - твоя функция из условия, x0 - это x^k-1, antigrad(x0) - это все из формулы (5.2) книга 
		x = x0 + kappa * antigrad(x0);
		x0 = x;
		++num_of_iter; //счетчик, num_of_calc - тоже не обращай внимания
		++num_of_grad;
		write.add_data(x0, foo(x0));
	}
	write.close();
	return { {x0.first,x0.second},foo(x0) };
}
std::pair<std::pair<double, double>, double> method_sopr_grad_for_kvad(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{

	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	std::pair<double, double> x = { };
	vec p(0, 0);
	double kappa  = mzs(eps/1000, foo, x0, antigrad(x0), num_of_calc);
	x = x0 + kappa * antigrad(x0);
	write.add_point(x);
	++num_of_grad;
	++num_of_iter;
	for (int i = 0; i < n-1; ++i)
	{
		p = (-dot_with_matrix_Q(antigrad(x0), antigrad(x)) / dot_with_matrix_Q(antigrad(x0), antigrad(x0))) *
			antigrad(x0) + antigrad(x);
		++num_of_grad;
		++num_of_iter;
		x0 = x;
		kappa = mzs(eps / 1000, foo, x0, p, num_of_calc);
		x = x0 + kappa * p;
		write.add_point(x);
	}

	write.close();
	return{ { x.first,x.second }, foo(x) };
}
std::pair<std::pair<double, double>, double> method_fletchera_rivsa(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	if (antigrad(x0).norm() < eps) return { x0,foo(x0) };
	++num_of_grad;
	++num_of_iter;
	std::pair<double, double> x = { };
	int k = 1;
	vec p(0, 0);
	p = antigrad(x0);
	double kappa;
	while (true)
	{
		kappa = mzs(eps / 1000, foo, x0, p, num_of_calc);
		x = x0 + kappa * p;
		write.add_point(x);
		if (antigrad(x).norm() < eps) break;
		if (k % 2 == 0)
		{
			p = antigrad(x);
			++k;
		}
		else 
		{
			p = (dot(antigrad(x) , antigrad(x)) / dot(antigrad(x0), p)) * p + antigrad(x);
			++k;
		}
		x0 = x;
		++num_of_grad;
		++num_of_iter;
	}
	write.close();
	return{ { x.first,x.second }, foo(x) };
}
std::pair<std::pair<double, double>, double> method_polaka_rivera(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	if (antigrad(x0).norm() < eps) return { x0,foo(x0) };
	std::pair<double, double> x = { };
	int k = 1;
	vec p(0, 0);
	p = antigrad(x0);
	++num_of_grad;
	++num_of_iter;
	double kappa;
	if (antigrad(x0).norm() < eps) return { x0, foo(x0) };
	while (true)
	{
		kappa = mzs(eps / 1000, foo, x0, p, num_of_calc);
		x = x0 + kappa * p;
		write.add_point(x);
		if (antigrad(x).norm() < eps) break;
		if (k % 2== 0)
		{
			p = antigrad(x);
			++k;
		}
		else 
		{
			p = (dot(antigrad(x) - antigrad(x0), antigrad(x)) / dot(antigrad(x0), antigrad(x0))) * antigrad(x0) + antigrad(x);
			++k;
		}
		x0 = x;
		++num_of_grad;
		++num_of_iter;
	}
	write.close();
	return{ { x.first,x.second }, foo(x) };
}
std::pair<std::pair<double, double>, double> method_sopr_grad_for_roz(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{

	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	if (antigrad(x0).norm() < eps) return { x0,foo(x0) };
	std::pair<double, double> x = { };
	int k = 1;
	vec p(0, 0);
	p = antigrad(x0);
	double kappa;
	++num_of_grad;
	++num_of_iter;
	if (antigrad(x0).norm() < eps) return { x0, foo(x0) };
	while (true)
	{
		kappa = mzs(eps / 1000, foo, x0, p, num_of_calc);
		x = x0 + kappa * p;
		write.add_point(x);
		if (antigrad(x).norm() < eps) break;
		if (k % 2 == 0)
		{
			p = antigrad(x);
		++k;
		}
		else 
		{
			p = (dot(antigrad(x), antigrad(x))/dot(antigrad(x0),p)) * p + antigrad(x);
			++k;
		}
		x0 = x;
		++num_of_grad;
		++num_of_iter;
	}
	write.close();
	return{ { x.first,x.second }, foo(x) };
}
std::pair<std::pair<double, double>, double> method_newtona(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>), matrix (*matr)(std::pair<double, double> point))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	write.add_level_line(foo(x0));
	std::pair<double, double> x;
	std::vector<double> pp(2);
	vec p(0, 0);
	matrix H = matr(x0);
	//{
	//	std::cout << "-----------" << std::endl;
	//	std::cout << "H=  " << std::endl;
	//	std::cout << H.get_A()[0][0] << "    " << H.get_A()[0][1] << std::endl;
	//	std::cout << H.get_A()[1][0] << "    " << H.get_A()[1][1] << std::endl;
	//	std::cout << "-----------" << std::endl;
	//}
	double lambda = lambda0;
	matrix E(base_type_of_matrix::singular, 2);
	while (antigrad(x0).norm() > eps)
	{
		//p = inverse(x0) * antigrad(x0) ;
		/*while (!H.is_positive())
		{
			H = H + E * lambda;
			lambda /= omega;
		}*/
		p = get_decision_gaus_(H, { antigrad(x0).get_x(), antigrad(x0).get_y() });
	//	p.set_x(pp[0]); p.set_y(pp[1]);
		++num_of_iter;
		++num_of_grad;
		num_of_second_d += 3;
		x = x0 + p;		
		x0 = x;
		H = matr(x0);
		write.add_point(x0);
		write.add_level_line(foo(x0));
//		lambda = lambda0;
		/*{
			std::cout << "-----------" << std::endl;
			std::cout << "H=  " << std::endl;
			std::cout << H.get_A()[0][0] << "    " << H.get_A()[0][1] << std::endl;
			std::cout << H.get_A()[1][0] << "    " << H.get_A()[1][1] << std::endl;
			std::cout << "-----------" << std::endl;
		}*/
	}
	write.close();
	return { {x0.first, x0.second},foo(x0) };
}
std::pair<std::pair<double, double>, double> method_newtona_naisk(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>), matrix(*matr)(std::pair<double, double> point))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	write.add_level_line(foo(x0));
	std::pair<double, double> x;
	vec p(0, 0);
	double kappa = kappa0;
	matrix H = matr(x0);
	double lambda = lambda0;
	matrix E(base_type_of_matrix::singular, 2);
	std::vector<double> pp(2);
	while (antigrad(x0).norm() > eps)
	{
	//	p = inverse(x0) * antigrad(x0);
	/*	while (!H.is_positive())
		{
			H = H + E * lambda;
			lambda /= omega;
		}*/
		p = get_decision_gaus_(H, { antigrad(x0).get_x(), antigrad(x0).get_y() });
		//p.set_x(pp[0]); p.set_y(pp[1]);
		kappa = mzs(eps / 1000, foo, x0, p, num_of_calc);
		++num_of_iter;
		++num_of_grad;
		num_of_second_d += 3;
		x = x0 + kappa * p;
		x0 = x;
		H = matr(x0);
//		lambda = lambda0*omega;
		write.add_point(x0);
		write.add_level_line(foo(x0));
	}
	write.close();
	return { {x0.first, x0.second},foo(x0) };
}
std::pair<std::pair<double, double>, double> method_markvarda(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>), matrix(*matr)(std::pair<double, double> point))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	write.add_level_line(foo(x0));
	std::pair<double, double> x;
	std::vector<double> pp(2);
	vec p(0, 0);
	double f1;
	double f2;
	double kappa = kappa0;
	double lambda = lambda0;
	matrix H = matr(x0);
	f1 = foo(x0);
	++num_of_calc;
	matrix E(base_type_of_matrix::singular,2);
	while (antigrad(x0).norm() > eps)
	{
		while (!H.is_positive())
		{
		//	std::cout << "+" << std::endl;
			H = H + E * lambda;
			lambda /= omega;
		}
		p = get_decision_gaus_(H, { antigrad(x0).get_x(), antigrad(x0).get_y() });
		//p.set_x(pp[0]); p.set_y(pp[1]);
		x = x0 + kappa*p;
		f2 = foo(x);
		while (f1 - f2 < omega * kappa * dot(antigrad(x0), p))
		{
			kappa = nu * kappa;
			x = x0 + kappa * p;
			f2 = foo(x);
			++num_of_calc;
		}
		++num_of_iter;
		++num_of_grad;
		++num_of_calc;
		num_of_second_d += 3;
		x0 = x;
		H = matr(x0);
		lambda = lambda0*omega;
		write.add_point(x0);
		write.add_level_line(foo(x0));
	//	lambda = lambda0 * omega;
	}
	return { {x0.first,x0.second}, foo(x0)};
}
std::pair<std::pair<double, double>, double> method_dfp(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	write.add_level_line(foo(x0));
	std::pair<double, double> x;
	vec p(0, 0);
	double kappa=kappa0;
	int k = 1;

	vec delt_x(0,0);
	vec delt_grad(0,0);
	matrix d_x(delt_x);
	matrix d_grad(delt_grad);

	matrix A(base_type_of_matrix::singular, 2);
	while (antigrad(x0).norm() > eps)
	{
		p = A * antigrad(x0);
		kappa = mzs(eps / 1000, foo, x0, p, num_of_calc);
		x = x0 + kappa * p;
		if (k % 2 == 0)
		{
			A = matrix(base_type_of_matrix::singular, 2);
			++k;
		}
		else
		{
		    delt_x = vec( x.first - x0.first , x.second - x0.second);
			delt_grad = vec(antigrad(x).get_x()- antigrad(x0).get_x(), antigrad(x).get_y()- antigrad(x0).get_y());
		    d_x = matrix(delt_x);
		    d_grad = matrix(delt_grad);
			A = A - ((d_x*d_x.transponate()) * (1./dot(delt_x,delt_grad))) - ((A* d_grad * d_grad.transponate()*A.transponate()) * (1./dot(A*delt_grad,delt_grad)));
			++k;
		}
		x0 = x;
		write.add_point(x0);
		write.add_level_line(foo(x0));
		++num_of_grad;
		++num_of_iter;
	}
	return  { {x0.first,x0.second}, foo(x0) };
}
std::pair<std::pair<double, double>, double> method_bfsh(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	write.add_level_line(foo(x0));
	std::pair<double, double> x;
	vec p(0, 0);
	double kappa = kappa0;
	int k = 1;

	vec delt_x(0, 0);
	vec delt_grad(0, 0);
	double rho = 0;
	vec r(0, 0);
	matrix d_x(delt_x);
	matrix d_grad(delt_grad);
	matrix R(r);

	matrix A(base_type_of_matrix::singular, 2);
	while (antigrad(x0).norm() > eps)
	{
	//	p = get_decision_gaus_(A, antigrad(x0));
		p = A * antigrad(x0);
		kappa = mzs(eps / 1000, foo, x0, p, num_of_calc);
		x = x0 + kappa * p;
		if (k % 2 == 0)
		{
			A = matrix(base_type_of_matrix::singular, 2);
			++k;
		}
		else
		{
			delt_x = vec(x.first - x0.first, x.second - x0.second);
			delt_grad = vec(antigrad(x).get_x() - antigrad(x0).get_x(), antigrad(x).get_y() - antigrad(x0).get_y());
			
			d_x = matrix(delt_x);
			d_grad = matrix(delt_grad);
			rho = dot(A * delt_grad, delt_grad);
			r = (A * delt_grad) * (1. / rho) - delt_x * (1. / dot(delt_x, delt_grad));
		
			R = matrix(r);
			A = A - ((d_x * d_x.transponate()) * (1. / dot(delt_x, delt_grad))) - ((A * d_grad * d_grad.transponate() * A.transponate()) * (1. / rho))+R*R.transponate()*rho;
			++k;
		}
		x0 = x;
		write.add_point(x0);
		write.add_level_line(foo(x0));
		++num_of_grad;
		++num_of_iter;
	}
	return  { {x0.first,x0.second}, foo(x0) };
}
std::pair<std::pair<double, double>, double> method_paulla(double(*foo)(std::pair<double, double>), vec(*antigrad)(std::pair<double, double>))
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	write.add_point(x0);
	write.add_level_line(foo(x0));
	std::pair<double, double> x;
	vec p(0, 0);
	double kappa = kappa0;
	int k = 1;

	vec delt_x(0, 0);
	vec delt_grad(0, 0);
	
	vec x_tild(0, 0);
	matrix d_x(delt_x);
	matrix d_grad(delt_grad);
	matrix R(x_tild);

	matrix A(base_type_of_matrix::singular, 2);
	while (antigrad(x0).norm() > eps)
	{
		//p = get_decision_gaus_(A, antigrad(x0));
		p = A* antigrad(x0);
		kappa = mzs(eps / 1000, foo, x0, p, num_of_calc);
		x = x0 + kappa * p;
		if (k % 2 == 0)
		{
			A = matrix(base_type_of_matrix::singular, 2);
			++k;
		}
		else
		{
			delt_x = vec(x.first - x0.first, x.second - x0.second);
			delt_grad = vec(antigrad(x).get_x() - antigrad(x0).get_x(), antigrad(x).get_y() - antigrad(x0).get_y());

			d_x = matrix(delt_x);
			d_grad = matrix(delt_grad);
			x_tild = delt_x + A * delt_grad;
			R = matrix(x_tild);

			matrix f1 = (R * R.transponate());
			double f2 = dot(x_tild, delt_grad);
			A = A - f1 / f2;
			++k;
		}
		x0 = x;
		write.add_point(x0);
		write.add_level_line(foo(x0));
		++num_of_grad;
		++num_of_iter;
	}
	return  { {x0.first,x0.second}, foo(x0) };
}
std::pair<std::pair<double, double>, double> method_pokord_spuska(double(*foo)(std::pair<double, double>))
{
//	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
//	write.add_point(x0);
//	write.add_level_line(foo(x0));
	std::pair<double, double> x;

	vec e1 = { 1, 0 };
	vec e2 = { 0, 1 };
	vec a = { 0,0 };

	while (true)
	{
		a.set_x(mzs(eps / 1000, foo, x0, e1, num_of_calc));

		x.first = x0.first + a.get_x();
		x.second = x0.second;
//		write.add_point(x);

		a.set_y(mzs(eps / 1000, foo, x, e2, num_of_calc));

		x.second = x0.second + a.get_y();
//		write.add_point(x);
	//	write.add_level_line(foo(x));

		++num_of_iter;
		if ((fabs(foo(x) - foo(x0)) < eps)&&(fabs(x-x0)<eps)) break;
		x0 = x;
	}
	x0 = x;
//	write.close();
	return  { {x0.first,x0.second}, foo(x0) };
}
std::pair<std::pair<double, double>, double> method_huka_jivsa(double(*foo)(std::pair<double, double>))
{
	//file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
 //   write.add_point(x0);
//	write.add_level_line(foo(x0));
	std::pair<double, double> x;
	vec e1 = { 1, 0 };
	vec e2 = { 0, 1 };
	vec b = { 1,1 };
	
	double left_f = 0;
	double right_f = 0;
	double top_f = 0;
	double down_f = 0;
	double mid_f = 0;
	std::map<double, std::pair<double, double> > line_f;
	std::map<double, std::pair<double, double>> row_f;
	double gamma = 2;
	vec p = { 0,0 };

	double a = 2;
	while (true)
	{
		
		num_of_calc += 5;
	
		line_f.insert({ foo(x0 - b.get_x() * e1), x0 - b.get_x() * e1 });
		line_f.insert({ foo(x0 + b.get_x() * e1), x0 + b.get_x() * e1 });
		line_f.insert({ foo(x0), x0  });
		std::pair<double, double> temp = { line_f.begin()->second.first, x0.second };
		
		row_f.insert({ foo(temp - b.get_y() * e2), temp - b.get_y() * e2 });
		row_f.insert({ foo(temp + b.get_y() * e2), temp + b.get_y() * e2 });
		row_f.insert({ foo(temp), temp }); 
		
		x = { line_f.begin()->second.first ,row_f.begin()->second.second };
		line_f.clear();
		row_f.clear();
		if (x == x0) { b = b * (1 / gamma); continue; }
		if ((x - x0 < eps) && (fabs(foo(x0) - foo(x)))) break;
		p = {x.first-x0.first, x.second-x0.second};
		a =  mzs(eps / 1000, foo, x, p, num_of_calc);
//	    write.add_point(temp);
//		write.add_point(x);
//		write.add_level_line(foo(x));
		x0 = x + a * p;
	//	++num_of_iter;
//		write.add_point(x0);
//		write.add_level_line(foo(x0));
	}
	x0 = x;
	//write.add_point(x0);
	//write.add_level_line(foo(x0));
	//write.close();
	return  { {x0.first,x0.second}, foo(x0) };

}
std::pair<std::pair<double, double>, double> method_rozenbroka(double(*foo)(std::pair<double, double>))
{
	/*file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));*/
	std::pair<double, double> x0 = nachalnaya_tochka;
//	write.add_point(x0);
//	write.add_level_line(foo(x0));
	std::pair<double, double> x=x0;
	vec p1 = { 1, 0 };
	vec p2 = { 0, 1 };
	vec a1 = { 0, 0 };
	vec a2 = { 0, 0 };
	vec b1 = { 0, 0 };
	vec b2 = { 0, 0 };
	double kappa1=0;
	double kappa2=0;
	while (true)
	{
		kappa1 = mzs(eps / 1000, foo, x, p1, num_of_calc);
		x = x + kappa1 * p1;
		kappa2= mzs(eps / 1000, foo, x, p2, num_of_calc);
		x = x + kappa2 * p2;
		if ((x - x0 < eps)&&(fabs(foo(x0)-foo(x)<eps))) break;
		if (kappa1 == 0) a1 = p1; else a1 = kappa1 * p1 + kappa2 * p2;
		if (kappa2 == 0) a2 = p2; else a2 = kappa2 * p2;
		b1 = a1;
		b2 = a2 - (dot(a2, b1)/ dot(b1, b1)) * b1;
		p1 = b1 * (1. / b1.norm());
		p2 = b2 * (1. / b2.norm());
		x0 = x;
	//	write.add_point(x0);
	//	write.add_level_line(foo(x0));
		++num_of_iter;
	}
//	write.close();
	return  { {x0.first,x0.second}, foo(x0) };
}
std::pair<std::pair<double, double>, double> method_regular_simplex(double(*foo)(std::pair<double, double>))
{
//	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
//	write.add_point(x0);
//	write.add_level_line(foo(x0));
	simplex simp(foo, x0);
//	write.add_simplex(simp);
	++num_of_calc;
	++num_of_calc;
	++num_of_calc;
	while (simp.max_length() > eps)
	{
		simp.next_step_for_regular(foo,num_of_calc);
//		write.add_simplex(simp);
		++num_of_iter;
		x0 = simp.trngl.begin()->second;
//		write.add_point(x0);
//		write.add_level_line(foo(x0));			
	}

//	write.close();
	return  { {x0.first,x0.second}, foo(x0) };
}
std::pair<std::pair<double, double>, double> method_irregular_simplex(double(*foo)(std::pair<double, double>))
{
//	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
//	write.add_point(x0);
//	write.add_level_line(foo(x0));
	simplex simp(foo, x0);
//	write.add_simplex(simp);	
	
	while (simp.max_length() > eps)
	{
		if (num_of_iter % 50==0) simp = simplex(foo, x0);
		simp.next_step_for_irregular(foo, num_of_calc);
//		write.add_simplex(simp);
		++num_of_iter;
		++num_of_calc;
		x0 = simp.trngl.begin()->second;
//		write.add_point(x0);
//		write.add_level_line(foo(x0));
	}

//	write.close();
	return  { {x0.first,x0.second}, foo(x0) };
}
std::pair<std::pair<double, double>, double> method_vnutr_shtraf(double(*foo)(std::pair<double, double>), int kvad_or_roz)
{
	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	std::pair<double, double> x = {};
	write.add_point(x0);
	write.add_level_line(foo(x0));
	if (kvad_or_roz == 1) {
		while (true)
		{

			x = method_huka_jivsa(vnutr_shtraf_for_kvad2).first;
			r_for_vnutr = r_for_vnutr / 1.4;
			if (fabs(foo(x0) - foo(x)) < eps) break;
			x0 = x;
			++num_of_iter;
			++num_of_calc;
			write.add_point(x0);
			write.add_level_line(vnutr_shtraf_for_kvad2(x0));
		}
	} else 
		while (true)
		{

			x = method_huka_jivsa(vnutr_shtraf_for_roz2).first;
			r_for_vnutr = r_for_vnutr / 1.4;
			if (fabs(foo(x0) - foo(x)) < eps) break;
			x0 = x;
			++num_of_iter;
			++num_of_calc;
			write.add_point(x0);
			write.add_level_line(vnutr_shtraf_for_roz2(x0));
		}
	r_for_vnutr = 32;
	write.close();
	return  { {x0.first,x0.second}, foo(x0) };
}
std::pair<std::pair<double, double>, double> method_vnesh_shtraf(double(*foo)(std::pair<double, double>), int kvad_or_roz)
{

	file_helper write(reinterpret_cast<void*>(foo), reinterpret_cast<void*>(kvad));
	std::pair<double, double> x0 = nachalnaya_tochka;
	std::pair<double, double> x = {};
	write.add_point(x0);
	write.add_level_line(foo(x0));
	if (kvad_or_roz == 1) {
		while (true)
		{

			x = method_huka_jivsa(vnesh_shtraf_for_kvad2).first;
			r_for_vnesh = r_for_vnesh / 0.5;
			if (fabs(foo(x0) - foo(x)) < eps) break;
			x0 = x;
			++num_of_iter;
			++num_of_calc;
			write.add_point(x0);
			write.add_level_line(vnesh_shtraf_for_kvad2(x0));
		}
 }
	else
		while (true)
		{

			x = method_huka_jivsa(vnesh_shtraf_for_roz2).first;
			r_for_vnesh = r_for_vnesh / 0.5;
			if (fabs(foo(x0) - foo(x)) < eps) break;
			x0 = x;
			++num_of_iter;
			++num_of_calc;
			write.add_point(x0);
			write.add_level_line(vnesh_shtraf_for_roz2(x0));
		}
	r_for_vnesh = 1;
	write.close();
	return  { {x0.first,x0.second}, foo(x0) };
}


double dot(const vec& antigrad1, const vec& antigrad2)
{

	return antigrad1.get_x() * antigrad2.get_x() + antigrad1.get_y() * antigrad2.get_y();
}
double dot_with_matrix_Q(const vec& antigrad1, const vec& antigrad2)
{

	vec new_vec(Q.get_A()[0][0] * antigrad1.get_x() + Q.get_A()[0][1] * antigrad1.get_y(),
		Q.get_A()[1][0] * antigrad1.get_x() + Q.get_A()[1][1] * antigrad1.get_y());
	return dot(new_vec, antigrad2);
}
double dot_with_matrix_H(const vec& antigrad1, const vec& antigrad2, std::pair<double, double> point)
{
	matrix H({ {12 * alpha * point.first * point.first - 4 * alpha * point.second + 2, -4 * alpha * point.first},{-4 * alpha * point.first, 2. * alpha} },2);
	vec new_vec(H.get_A()[0][0] * antigrad1.get_x() + H.get_A()[0][1] * antigrad1.get_y(),
		H.get_A()[1][0] * antigrad1.get_x() + H.get_A()[1][1] * antigrad1.get_y());
	return dot(new_vec, antigrad2);
}
matrix get_H_for_roz(std::pair<double, double> point)
{
	matrix H({ {2 + 12 * alpha * point.first * point.first - 4 * alpha * point.second,
		-4 * alpha * point.first},{-4 * alpha * point.first, 2. * alpha} },2);
	return H;
}
matrix get_H_for_roz_inverse(std::pair<double, double> point)
{
	matrix H({ {1 / (2 + 4 * alpha * (point.first * point.first - point.second)),

		point.first / (1 + 2 * alpha * (point.first * point.first - point.second))},

		{point.first / (1 + 2 * alpha * (point.first * point.first - point.second)),

		1. / (2 * alpha) + (2 * point.first * point.first) /
		(1 + 2 * alpha * point.first * point.first - 2 * alpha * point.second)} },2);
	return H;
}
matrix get_H_for_kvad(std::pair<double, double> point)
{
	matrix H({ {4,-4},{-4,10} },2);
	return H;
}
matrix get_H_for_kvad_inverse(std::pair<double, double> point)
{
	matrix H({ {5. / 12,1. / 6},{1. / 6,1. / 6} },2);
	return H;
}
double operator-(std::pair<double, double> p1, std::pair<double, double> p2)
{
	std::pair<double, double> p3 = { p1.first - p2.first,p1.second - p2.second };
	double res = sqrt(p3.first * p3.first + p3.second * p3.second);
	return res;
}

vec operator*(matrix A, vec b)
{
	std::vector<double> res(n);
	for (int i = 0; i < n; ++i)
	{
		res[i] += A.get_A()[i][0] * b.get_x();
		res[i] += A.get_A()[i][1] * b.get_y();
	}
	vec answ(res[0], res[1]);
	return answ;
}

double  vnutr_shtraf_for_kvad1(std::pair<double, double> x)
{
	if ((x.first > 0) && (x.second > 2) && (x.second + x.first < 10))
		return (kvad(x) - r_for_vnutr * (1 / (-x.first) + 1 / (2-x.second) + 1 / (x.first + x.second - 10)));
	else return (kvad(x) + DBL_MAX);
}
double  vnutr_shtraf_for_roz1(std::pair<double, double> x)
{
	if ((x.first > 0) && (x.second > 2) && (x.second + x.first < 10))
		return (roz(x) - r_for_vnutr * (1 / (-x.first) + 1 / (2-x.second) + 1 / (x.first + x.second - 10)));
	else return (roz(x) + DBL_MAX);
}
double  vnutr_shtraf_for_kvad2(std::pair<double, double> x)
{
	if ((((x.first + 3) * (x.first +3)) / 4 + ((x.second + 4) * (x.second + 4)) / 9 - 5) < 0)
		return (kvad(x) - r_for_vnutr * 1 / (((x.first + 3) * (x.first + 3)) / 4 + ((x.second + 4) * (x.second + 4)) / 9 - 5));
	else return (kvad(x) + DBL_MAX);
}
double  vnutr_shtraf_for_roz2(std::pair<double, double> x)
{
	if ((((x.first + 3) * (x.first + 3)) / 4 + ((x.second + 4) * (x.second + 4)) / 9 - 5) < 0)
	return (roz(x) - r_for_vnutr * 1 / (((x.first + 3) * (x.first + 3)) / 4 + ((x.second + 4) * (x.second + 4)) / 9 - 5));
	else return (roz(x) + DBL_MAX);
}

double  vnesh_shtraf_for_kvad1(std::pair<double, double> x)
{
	if ((x.first > 0) && (x.second > 2) && (x.second + x.first < 10))
		return kvad(x);
	else
	{
		double temp_x = -x.first;
		double temp_y = 2-x.second;
		double temp_x_y = x.second + x.first-10;
		if (temp_x <= 0) temp_x = 0;
		if (temp_y <= 0) temp_y = 0;
		if (temp_x_y <= 0) temp_x_y = 0;
		return (kvad(x) + r_for_vnesh * (temp_x*temp_x + temp_y*temp_y + temp_x_y* temp_x_y));
	}		
}
double  vnesh_shtraf_for_roz1(std::pair<double, double> x)
{
	if ((x.first > 0) && (x.second > 2) && (x.second + x.first < 10))
		return roz(x);
	else
	{
		double temp_x = -x.first;
		double temp_y = 2-x.second;
		double temp_x_y = x.second + x.first-10;
		if (temp_x <= 0) temp_x = 0;
		if (temp_y <= 0) temp_y = 0;
		if (temp_x_y <= 0) temp_x_y = 0;
		return (roz(x) + r_for_vnesh * (temp_x * temp_x + temp_y * temp_y + temp_x_y * temp_x_y));
	}
}
double  vnesh_shtraf_for_kvad2(std::pair<double, double> x)
{
	if ((((x.first + 3) * (x.first + 3)) / 4 + ((x.second + 4) * (x.second + 4)) / 9 - 5) < 0)
		return kvad(x);
	else
	{
		double temp_x = (((x.first + 3) * (x.first + 3)) / 4 + ((x.second + 4) * (x.second + 4)) / 9 - 5);
		if (temp_x <= 0) temp_x = 0;
		return (kvad(x) + r_for_vnesh * temp_x * temp_x);
	}
}
double  vnesh_shtraf_for_roz2(std::pair<double, double> x)
{
	if ((((x.first + 3) * (x.first + 3)) / 4 + ((x.second + 4) * (x.second + 4)) / 9 - 5) < 0)
		return roz(x);
	else
	{
		double temp_x = (((x.first + 3) * (x.first + 3)) / 4 + ((x.second + 4) * (x.second + 4)) / 9 - 5);
		if (temp_x <= 0) temp_x = 0;
		return (roz(x) + r_for_vnesh * temp_x * temp_x);
	}
}