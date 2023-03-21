
#include <format>
#include <iostream>
#include <cmath>

using func = double(*)(double);

const double h = 0.01;
double sign(double x)
{
	if (x != x)
		return x;
	return std::copysign(1, x) * (x != 0);
}
void report_approximation(size_t step, double x, double y)
{
	std::cout << std::format("x{} = {:<+22.16g} y{} = {:<+22.16g}\n", step, x, step, y);
}
double get_second_difference(func f, double x)
{
	return f(x + h) - 2 * f(x) + f(x - h);
}
bool sign_matches(double a, double b)
{
	return std::copysign(a, b) == a;
}
double finite_difference_derivative(func f, double x)
{
	return (f(x + h) - f(x - h)) / (2 * h);
}
double finite_difference_second_derivative(func f, double x)
{
	return get_second_difference(f, x) / (2 * h);
}


double run_secant_method(func f, double x_precision, double x1, double x2, size_t max_steps)
{
	size_t step = 0;
	while (true)
	{
		if (step == max_steps)
			return NAN;

		const double f1 = f(x1);
		const double f2 = f(x2);
		const double dx = f(x1) * (x2 - x1) / (f2 - f1);
		
		const double x3 = x1 - dx;
		report_approximation(++step, x3, f(x3));

		if (std::isnan(x3))
			return NAN;

		if (std::abs(dx) <= x_precision / 2)
			return x3;

		x1 = x2;
		x2 = x3;
	}
}
double run_chord_method(func f, double x_precision, double x1, double x2, size_t max_steps)
{
	size_t step = 0;
	if (!sign_matches(f(x1), get_second_difference(f, x1)))
		std::swap(x1, x2);

	while (true)
	{
		if (step == max_steps)
			return NAN;

		const double f1 = f(x1);
		const double f2 = f(x2);
		const double dx = f(x2) * (x2 - x1) / (f2 - f1);

		const double x3 = x2 - dx;
		report_approximation(++step, x3, f(x3));

		if (std::isnan(x3))
			return NAN;

		if (std::abs(dx) <= x_precision / 2)
			return x3;

		x2 = x3;
	}
}
double run_dichotomy_method(func f, double x_precision, double x1, double x2, size_t max_steps)
{
	const double s1 = sign(f(x1));
	const double s2 = sign(f(x2));
	if (s1 == 0)
		return x1;
	if (s2 == 0)
		return x2;
	if (s1 == s2)
		return NAN;

	if (s1 != 1)
		std::swap(x1, x2);

	size_t step = 0;
	while (true)
	{
		if (step++ == max_steps)
			return NAN;

		const double mid = (x1 + x2) / 2;

		const double length = std::abs(x2 - x1);
		if (length <= x_precision)
			return mid;

		const double mid_val = f(mid);
		report_approximation(step, mid, mid_val);

		const double mid_sign = sign(mid_val);

		if (mid_sign == 0)
			return mid;

		if (mid_sign == 1)
			x1 = mid;
		else if (mid_sign == -1)
			x2 = mid;
		else
			return NAN;
	}
}
double run_newthon_method(func f, double x_precision, double x, size_t max_steps)
{
	size_t step = 0;
	while (true)
	{
		if (step++ == max_steps)
			return NAN;
		if (std::isnan(x))
			return NAN;

		const double dx = f(x) / finite_difference_derivative(f, x);
		x = x - dx;
		report_approximation(step, x, f(x));


		if (std::abs(dx) <= x_precision / 2)
			return x;
		if (std::isnan(x))
			return NAN;
	}
}
double run_halley_method(func f, double x_precision, double x, size_t max_steps)
{
	size_t step = 0;
	while (true)
	{
		if (step++ == max_steps)
			return NAN;
		if (std::isnan(x))
			return NAN;

		const double dfdx = finite_difference_derivative(f, x);
		const double a = f(x) / dfdx;
		const double b = (1 - a * finite_difference_second_derivative(f, x) / (2 * dfdx));
		const double dx = a / b;
		x = x - dx;
		report_approximation(step, x, f(x));


		if (std::abs(dx) <= x_precision / 2)
			return x;
		if (std::isnan(x))
			return NAN;
	}
}

double f(double x)
{
	return std::sin(x) / x + x;
}

int main()
{
	double x = run_halley_method(f, 0, 1, 100);
}
