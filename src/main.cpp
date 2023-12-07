//#include <iostream>
#include <iomanip>
#include "necessary_func.h"
#include "simplex.h"


using namespace std;

extern const int alpha =4;
extern const int n = 2;
int num_of_iter = 0;
int num_of_calc = 0;
int num_of_grad = 0;
int num_of_second_d = 0;
double r_for_vnutr = 32;
double r_for_vnesh = 1;
std::string name_of_current_foo="";
extern  const  double kappa0 = 1.;
extern const double omega = 0.5;
extern const double lambda0 = 1.;
extern const double nu = 0.5;
extern const double eps = 0.01;
extern  const std::pair<double, double> nachalnaya_tochka = {-2,-8};

extern const matrix Q = matrix({ {1,-1}, {-1, 2.5} }, 2);



void choose_metod(int type_of_method_);
void write(pair < pair<double, double>, double> res);
int main()
{

	int type_of_method =18;
	// 1 - лернд дпнакемхъ ьюцю
	// 2 - лернд мюхяйнпеиьецн яосяйю
	// 3 - лернд янопъфеммшу цпюдхемрнб
	// 4 - лернд ткервепю-пхбяю
	// 5 - лернд онкюйю-пхбэепю
	// 6 - лернд мэчрнмю
	// 7 - лернд мэчрнмю C мюхяйнпеиьхл яосяйнл
	// 8 - лернд люпйбюпдрю
	// 9 - лернд дто
	// 10 - лернд ать
	// 11 - лернд оюсщккю
	// 12 - лернд жхйкхвеяйнцн онйннпдхмюрмнцн яосяйю
	// 13 - лернд усйю-дфхбяю
	// 14 - лернд пнгемапнйю
	// 15 - лернд пецскъпмнцн яхлокейяю
	// 16 - лернд мепецскъпмнцн яхлокейяю
	// 17 - лернд бмсрпеммху ьрпютмшу тсмйжхи
	// 18 - лернд бмсрпеммху ьрпютмшу тсмйжхи
	choose_metod(type_of_method);
}

void choose_metod(int type_of_method_)
{
	double(*foo_roz)(std::pair<double, double>) = roz;
	vec(*antigrad_roz)(std::pair<double, double>) = antigrad_for_roz;

	double(*foo_kvad)(std::pair<double, double>) = kvad;
	vec(*antigrad_kvad)(std::pair<double, double>) = antigrad_for_kvad;

	switch (type_of_method_)
	{
	case 1:
		cout << "KVAD FOO" << "\n" << "-------------------------"<<endl;
		write( method_droblenia_shaga(foo_kvad, antigrad_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_droblenia_shaga(foo_roz, antigrad_roz));
		break;
	case 2:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_naiscor_spuska (foo_kvad, antigrad_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_naiscor_spuska(foo_roz, antigrad_roz));
		break;
	case 3:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_sopr_grad_for_kvad (foo_kvad, antigrad_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_sopr_grad_for_roz (foo_roz, antigrad_roz));
		break;
	case 4:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_fletchera_rivsa (foo_kvad, antigrad_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_fletchera_rivsa(foo_roz, antigrad_roz));
		break;
	case 5:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_polaka_rivera (foo_kvad, antigrad_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_polaka_rivera(foo_roz, antigrad_roz));
		break;
	case 6:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		//write(method_newtona(foo_kvad, antigrad_kvad, get_H_for_kvad_inverse));
		write(method_newtona(foo_kvad, antigrad_kvad, get_H_for_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		//write(method_newtona(foo_roz, antigrad_roz, get_H_for_roz_inverse));
		write(method_newtona(foo_roz, antigrad_roz, get_H_for_roz));
		break;
	case 7:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		//write(method_newtona_naisk(foo_kvad, antigrad_kvad, get_H_for_kvad_inverse));
		write(method_newtona_naisk(foo_kvad, antigrad_kvad, get_H_for_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		//write(method_newtona_naisk(foo_roz, antigrad_roz, get_H_for_roz_inverse));
		write(method_newtona_naisk(foo_roz, antigrad_roz, get_H_for_roz));
		break;
	case 8:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_markvarda(foo_kvad, antigrad_kvad, get_H_for_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_markvarda(foo_roz, antigrad_roz, get_H_for_roz));
		break;
	case 9:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_dfp(foo_kvad, antigrad_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_dfp(foo_roz, antigrad_roz));
		break;
	case 10:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_bfsh(foo_kvad, antigrad_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_bfsh(foo_roz, antigrad_roz));
		break;
	case 11:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_paulla(foo_kvad, antigrad_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_paulla(foo_roz, antigrad_roz));
		break;
	case 12:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_pokord_spuska(foo_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_pokord_spuska(foo_roz));
		break;
	case 13:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_huka_jivsa(foo_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_huka_jivsa(foo_roz));
		break;
	case 14:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_rozenbroka(foo_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_rozenbroka(foo_roz));
		break;
	case 15:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_regular_simplex(foo_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_regular_simplex(foo_roz));
		break;
	case 16:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_irregular_simplex(foo_kvad));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_irregular_simplex(foo_roz));
		break;
	case 17:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_vnutr_shtraf(foo_kvad,1));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_vnutr_shtraf(foo_roz,2));
		break;
	case 18:
		cout << "KVAD FOO" << "\n" << "-------------------------" << endl;
		write(method_vnesh_shtraf(foo_kvad, 1));
		cout << "ROZENBROK FOO" << "\n" << "-------------------------" << endl;
		write(method_vnesh_shtraf(foo_roz, 2));
		break;
	default:
		break;
	}
}

void write(pair < pair<double, double>, double> res)
{
	cout << "num_of_iter = " << num_of_iter << endl;
	cout << "num_of_calc = " << num_of_calc << endl;
	cout << "num_of_grad = " << num_of_grad << endl;
	cout << "num_of_second_d = " << num_of_second_d << endl;
	cout << "x_min = " << setprecision(8) << "( " << res.first.first << " , " << res.first.second << " )" << endl;
	cout << "y_min = " << setprecision(8) << res.second << endl;
	num_of_iter = 1;
	num_of_calc = 1;
	num_of_grad = 1;
}