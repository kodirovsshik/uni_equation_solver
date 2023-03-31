
#include <format>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <ranges>


class formula
{
	std::string m_expr;

public:
	formula(const std::string& expression);
	bool validate(std::string& err_msg, const char*& err_pos) const noexcept;
	double operator()(double argument) const;
};
using func = const formula&;


const double h = 0.01;
double sign(double x)
{
	if (x != x)
		return x;
	return std::copysign(1, x) * (x != 0);
}
void report_approximation(size_t step, double x, double y)
{
	std::cout << std::format("x{:02} = {:<+22.16g} y{:02} = {:<+22.16g}\n", step, x, step, y);
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
bool isinfnan(double x)
{
	return std::isnan(x) || std::isinf(x);
}
double estimate_root(func f, double x1, double x2)
{
	const double mid = (x1 + x2) / 2;
	return std::ranges::min({ x1, x2, mid }, {}, [&](double x) { return std::abs(f(x)); });
}

double run_secant_method(func f, double x_precision, double x1, double x2, size_t max_steps = SIZE_MAX)
{
	std::cout << "\nSecant method:\n";
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

		if (isinfnan(x3))
			return NAN;

		if (std::abs(dx) <= x_precision / 2)
			return x3;

		x1 = x2;
		x2 = x3;
	}
}
double run_chord_method(func f, double x_precision, double x1, double x2, size_t max_steps = SIZE_MAX)
{
	std::cout << "\nChord method:\n";
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

		if (isinfnan(x3))
			return NAN;

		if (std::abs(dx) <= x_precision / 2)
			return x3;

		x2 = x3;
	}
}
double run_dichotomy_method(func f, double x_precision, double x1, double x2, size_t max_steps = SIZE_MAX)
{
	std::cout << "\nDichotomy method:\n";
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
double run_newthon_method(func f, double x_precision, double x, size_t max_steps = SIZE_MAX)
{
	std::cout << "\nNewthon method:\n";
	size_t step = 0;
	while (true)
	{
		if (step++ == max_steps)
			return NAN;
		if (isinfnan(x))
			return NAN;

		const double dx = f(x) / finite_difference_derivative(f, x);
		x = x - dx;
		report_approximation(step, x, f(x));


		if (std::abs(dx) <= x_precision / 2)
			return x;
		if (isinfnan(x))
			return NAN;
	}
}
double run_halley_method(func f, double x_precision, double x, size_t max_steps = SIZE_MAX)
{
	std::cout << "\nHalley method:\n";
	size_t step = 0;
	while (true)
	{
		if (step++ == max_steps)
			return NAN;
		if (isinfnan(x))
			return NAN;

		const double dfdx = finite_difference_derivative(f, x);
		const double a = f(x) / dfdx;
		const double b = (1 - a * finite_difference_second_derivative(f, x) / (2 * dfdx));
		const double dx = a / b;
		x = x - dx;
		report_approximation(step, x, f(x));


		if (std::abs(dx) <= x_precision / 2)
			return x;
		if (isinfnan(x))
			return NAN;
	}
}
double run_simple_iterations_method(func f, double x_precision, double x1, double x2, size_t max_steps = SIZE_MAX)
{
	std::cout << "\nSimple iterations method:\n";
	double lambda, x;
	{
		double dfdx1 = finite_difference_derivative(f, x1);
		double dfdx2 = finite_difference_derivative(f, x2);

		if (sign(dfdx1) != sign(dfdx2))
		{
			std::cout << "function rejected: derivative's sign alternates\n";
			return NAN;
		}

		double adfdx1 = std::copysign(dfdx1, 1);
		double adfdx2 = std::copysign(dfdx2, 1);

		const double max_derivative = std::max(adfdx1, adfdx2);
		if (max_derivative == 0)
			return NAN;
		lambda = std::copysign(1, dfdx1) / max_derivative;
		
		x = estimate_root(f, x1, x2);
		report_approximation(0, x, f(x));
	}

	size_t step = 0;
	while (true)
	{
		if (step++ == max_steps)
			return NAN;
		if (isinfnan(x))
			return NAN;

		const double y = f(x), dx = lambda * y;
		x -= dx;
		report_approximation(step, x, y);

		if (std::abs(dx) < x_precision / 2)
			return x;
	}
}


int main()
{
	formula f("");

	while (true)
	{
		std::cout << "f(x) = ";
		std::string str;
		std::getline(std::cin, str);

		f = formula(str);
		const char* perr = nullptr;
		if (f.validate(str, perr))
			break;
		std::cout << "error: " << str << "\nstarting at " << std::string_view(perr, std::min<size_t>(10, strlen(perr))) << "\n";
	}

	const double prec = 1e-8, l = -2, r = 3, x0 = (l + r) / 2;
	run_secant_method(f, prec, l, r, 100);
	run_chord_method(f, prec, l, r, 100);
	run_dichotomy_method(f, prec, l, r, 100);
	run_newthon_method(f, prec, x0, 100);
	run_halley_method(f, prec, x0, 100);
	run_simple_iterations_method(f, prec, l, r, 100);
}



#include <charconv>
#include <unordered_map>

struct parse_context
{
	parse_context(const std::string_view& str) noexcept;

	const char* p = nullptr;
	const char* pe = nullptr;
	double arg = NAN;
	std::string error_message;
	bool dry = false;
};

bool try_parse_literal(parse_context& cntxt, double& var);
bool try_parse_variable(parse_context& cntxt, double& var);
bool try_parse_function(parse_context& cntxt, double& var);
bool parse_expression(parse_context& cntxt, double& var);

template<class T>
bool try_match_dictionary(
	const std::unordered_map<std::string_view, T>& map, 
	parse_context& cntxt, T& result, std::string_view postfix = "")
{
	const size_t leftover_len = cntxt.pe - cntxt.p;
	for (auto&& [key, value] : map)
	{
		const size_t key_len = key.size();
		const size_t post_len = postfix.size();

		if (leftover_len < key_len + post_len)
			continue;

		if (std::string_view(cntxt.p, key_len) != key)
			continue;

		if (std::string_view(cntxt.p + key_len, post_len) != postfix)
			continue;

		cntxt.p += key_len + post_len;
		result = value;
		return true;
	}

	return false;
}

bool try_parse_literal(parse_context& cntxt, double& var)
{
	auto& p = cntxt.p;
	auto& pe = cntxt.pe;
	auto result = std::from_chars(p, pe, var);
	p = result.ptr;
	return result.ec == std::errc{};
}
bool try_parse_variable(parse_context& cntxt, double& var)
{
	auto& p = cntxt.p;
	auto& pe = cntxt.pe;
	auto& arg = cntxt.arg;
	if (p == pe)
		return false;
	const bool is_match = *p == 'x';
	if (is_match)
		var = arg;
	p += is_match;
	return is_match;
}
bool try_parse_function(parse_context& cntxt, double& var)
{
	using unary_function = double(*)(double);

	static const  std::unordered_map<std::string_view, unary_function> functions = {
		{ "sin", &sin },
		{ "cos", &cos },
		{ "tan", &tan },
		{ "ctg", [](double x) { return 1 / tan(x); }},
		{ "sqrt", &sqrt },
		{ "cbrt", &cbrt },
		{ "sqr", [](double x) { return x * x; }},
		{ "abs", &abs },
		{ "exp", &exp },
		{ "ln", &log },
		{ "lg", &log10 },
		{ "log2", &log2 },
	};

	auto& p = cntxt.p;
	auto& pe = cntxt.pe;
	const size_t leftover_length = pe - p;

	unary_function parsed_func_ptr;
	if (!try_match_dictionary(functions, cntxt, parsed_func_ptr, "("))
		return false;

	double func_arg;
	if (!parse_expression(cntxt, func_arg))
		return false;

	if (p == pe || *p != ')')
		return false;
	++p;

	if (cntxt.dry)
		var = 0;
	else
		var = parsed_func_ptr(func_arg);
	return true;
}
bool is_unary_operator(char c)
{
	return std::string_view("+-").find(c) != std::string::npos;
}
bool try_parse_unary_expr(parse_context& cntxt, double& var)
{
	auto& p = cntxt.p;
	auto& pe = cntxt.pe;
	
	if (p == pe)
		return false;

	char op = *p;
	if (!is_unary_operator(op))
		return false;
	
	++p;
	bool parse_ok = parse_expression(cntxt, var);
	if (parse_ok && op == '-')
		var = -var;
	return parse_ok;
}
bool try_parse_nested_expression(parse_context& cntxt, double& var)
{
	if (cntxt.p == cntxt.pe || *cntxt.p != '(')
		return false;

	++cntxt.p;
	if (!parse_expression(cntxt, var))
		return false;
	if (*cntxt.p != ')')
		return false;
	++cntxt.p;
	return true;
}
bool try_parse_operand(parse_context& cntxt, double& var)
{
	parse_context cntxt_copy = cntxt;

	if (try_parse_literal(cntxt, var))
		return true;
	cntxt = cntxt_copy;

	if (try_parse_variable(cntxt, var))
		return true;
	cntxt = cntxt_copy;

	if (try_parse_function(cntxt, var))
		return true;
	cntxt = cntxt_copy;
	
	if (try_parse_nested_expression(cntxt, var))
		return true;
	//cntxt = cntxt_copy;

	return false;
}

using binary_func = double(*)(double, double);
struct operand
{
	double value;
	binary_func next_operator = nullptr;
};

double op_plus(double a, double b) { return a + b; }
double op_minus(double a, double b) { return a - b; }
double op_multiply(double a, double b) { return a * b; }
double op_divide(double a, double b) { return a / b; }
double op_divide_integer(double a, double b) { return trunc(a / b * (1 + 2 * DBL_EPSILON)); }
double op_remainder(double a, double b) { return fmod(a, b); }
double op_power(double a, double b) { return pow(a, b); }

bool takes_precedence(binary_func op1, binary_func op2)
{
	if (op1 == op_power)
		return true;
	
	static const std::unordered_map<binary_func, int> operator_precedence = {
		{ op_plus, 1 },
		{ op_minus, 1 },
		{ op_multiply, 2 },
		{ op_divide, 2 },
		{ op_divide_integer, 2 },
		{ op_remainder, 2 },
		{ op_power, 3 },
	};

	return operator_precedence.at(op1) > operator_precedence.at(op2);
}

bool try_parse_binary_expr(parse_context& cntxt, double& var, operand* p_operand1 = nullptr)
{
	static const std::unordered_map<std::string_view, binary_func> binary_operators = {
		{"+", op_plus},
		{"-", op_minus},
		{"*", op_multiply},
		{"/", op_divide},
		{"//", op_divide_integer},
		{"%", op_remainder},
		{"^", op_power},
	};

	auto& p = cntxt.p;
	auto& pe = cntxt.pe;

	operand _operand1_local;

	if (p_operand1 == nullptr)
	{
		if (!try_parse_operand(cntxt, _operand1_local.value))
			return false;

		if (p == pe || *p == ')')
		{
			var = _operand1_local.value;
			return true;
		}

		if (!try_match_dictionary(binary_operators, cntxt, _operand1_local.next_operator))
			return false;

		p_operand1 = &_operand1_local;
	}

	operand& operand1 = *p_operand1;
	operand operand2;

	if (!try_parse_operand(cntxt, operand2.value))
		return false;
	if (p == pe || *p == ')')
	{
		if (cntxt.dry)
			var = 0;
		else
			var = operand1.next_operator(operand1.value, operand2.value);
		return true;
	}

	if (!try_match_dictionary(binary_operators, cntxt, operand2.next_operator))
		return false;

	if (takes_precedence(operand2.next_operator, operand1.next_operator))
	{
		double rest;
		if (!try_parse_binary_expr(cntxt, rest, &operand2))
			return false;
		var = operand1.next_operator(operand1.value, rest);
		return true;
	}
	else
	{
		operand2.value = operand1.next_operator(operand1.value, operand2.value);
		return try_parse_binary_expr(cntxt, var, &operand2);
	}
}

bool parse_expression(parse_context& cntxt, double& var)
{
	if (cntxt.p == cntxt.pe)
		return false;

	parse_context cntxt_copy = cntxt;

	if (try_parse_binary_expr(cntxt, var))
		return true;
	cntxt = cntxt_copy;

	if (try_parse_function(cntxt, var))
		return true;
	cntxt = cntxt_copy;

	if (try_parse_unary_expr(cntxt, var))
		return true;
	//cntxt = cntxt_copy;

	return false;
}

formula::formula(const std::string& expression)
{
	this->m_expr.reserve(expression.size());
	for (auto&& c : expression)
	{
		if (!std::isblank(c))
			this->m_expr += c;
	}
}

bool formula::validate(std::string& err_msg, const char*& err_pos) const noexcept
{
	int parentheses_level = 0;
	for (auto&& c : this->m_expr)
	{
		if (c == '(')
			++parentheses_level;
		if (c == ')')
			--parentheses_level;
		if (parentheses_level < 0)
		{
			err_pos = &c;
			err_msg = "unexpected ')'";
			return false;
		}
	}
	if (parentheses_level != 0)
	{
		err_pos = this->m_expr.data() + this->m_expr.size();
		err_msg = "no ')' to match '('";
		return false;
	}

	parse_context cntxt(this->m_expr);
	cntxt.dry = true;

	double unused;
	bool dry_parse_result = parse_expression(cntxt, unused);
	
	err_msg = std::move(cntxt.error_message);
	if (!dry_parse_result && err_msg.empty())
		err_msg = "failed to classify token sequence";
	err_pos = cntxt.p;
	return dry_parse_result;
}


double formula::operator()(double x) const
{
	double result = std::bit_cast<double>(0xFFFFFFFFFFFFFFFF);

	parse_context cntxt(this->m_expr);
	cntxt.arg = x;

	parse_expression(cntxt, result);
	return result;
}

parse_context::parse_context(const std::string_view& str) noexcept
	: p(str.data()), pe(str.data() + str.size())
{
}
