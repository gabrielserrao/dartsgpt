#include <complex>
#include "dartsflash/global/global.hpp"
#include "dartsflash/maths/geometry.hpp"
#include "dartsflash/maths/maths.hpp"
#include "dartsflash/maths/modifiedcholeskys99.hpp"
#include "dartsflash/maths/linesearch.hpp"
#include "dartsflash/maths/root_finding.hpp"

int test_combinations();
int test_simplex();
int test_function_pass();
int test_rootfinding();
int test_line_search();
int test_cubic_roots();

int main() 
{
    int error_output = 0;

    error_output += test_combinations();
    error_output += test_simplex();
    error_output += test_function_pass();
    error_output += test_rootfinding();
    error_output += test_line_search();
    error_output += test_cubic_roots();

    return error_output;
}

int test_combinations()
{
    int error_output = 0;

    std::vector<std::pair<Combinations, std::vector<std::vector<int>>>> combinations = {
        {Combinations(2, 2), {{0, 1}}},
        {Combinations(3, 2), {{0, 1}, {0, 2}, {1, 2}}},
        {Combinations(4, 2), {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}}},
        {Combinations(3, 3), {{0, 1, 2}}},
        {Combinations(4, 3), {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}},
        {Combinations(5, 3), {{0, 1, 2}, {0, 1, 3}, {0, 1, 4}, {0, 2, 3}, {0, 2, 4}, {0, 3, 4}, {1, 2, 3}, {1, 2, 4}, {1, 3, 4}, {2, 3, 4}}},
        {Combinations(4, 4), {{0, 1, 2, 3}}},
        {Combinations(5, 4), {{0, 1, 2, 3}, {0, 1, 2, 4}, {0, 1, 3, 4}, {0, 2, 3, 4}, {1, 2, 3, 4}}},
    };

    for (auto it: combinations)
    {
        for (size_t i = 0; i < it.first.combinations.size(); i++)
        {
            std::vector<int> combination = it.first.combinations[i];
            for (size_t ii = 0; ii < combination.size(); ii++)
            {
                if (combination[ii] != it.second[i][ii])
                {
                    error_output++;
                    std::cout << "Different combinations\n";
        			print("n_elements, combination_length", std::vector<int>{it.first.n_elements, it.first.combination_length});
        			print("combinations", it.first.combinations);
		        	print("ref", it.second);
                    break;
                }
            }
        }
    }

    if (error_output > 0)
	{
		print("Errors occurred in test_combinations()", error_output);
	}
	else
	{
		print("No errors occurred in test_combinations()", error_output);
	}
	return error_output;
}

int test_simplex()
{
    int error_output = 0;

    // Is between two points? (2D simplex/line)
    std::vector<std::vector<double>> coords = {{0.}, {1.}};
    std::vector<std::vector<double>> points = {{0.2}, {1.}, {1.2}};
    std::vector<bool> in_simplex = {true, true, false};
    for (size_t i = 0; i < points.size(); i++)
    {
        if (is_in_simplex(points[i], coords) != in_simplex[i])
        {
            error_output++;
            std::cout << (!in_simplex[i] ? "Point is in 2D simplex (line), but shouldn't be" : "Point isn't in 2D simplex (line), but should be") << "\n";
            print("Point", points[i]);
        }
    }

    coords = {{0.5}, {2.}};
    points = {{0.2}, {0.5}, {1.}, {2.}, {3.}};
    in_simplex = {false, true, true, true, false};
    for (size_t i = 0; i < points.size(); i++)
    {
        if (is_in_simplex(points[i], coords) != in_simplex[i])
        {
            error_output++;
            std::cout << (!in_simplex[i] ? "Point is in 2D simplex (line), but shouldn't be" : "Point isn't in 2D simplex (line), but should be") << "\n";
            print("Point", points[i]);
        }
    }

    // Is in triangle? (3D simplex/triangle)
    coords = {{0., 0.}, {0., 1.}, {1., 0.}};
    points = {{0.2, 0.2}, {0., 1.}, {0.5, 0.5}, {1., 1.}};
    in_simplex = {true, true, true, false};
    for (size_t i = 0; i < points.size(); i++)
    {
        if (is_in_simplex(points[i], coords) != in_simplex[i])
        {
            error_output++;
            std::cout << (!in_simplex[i] ? "Point is in 3D simplex (triangle), but shouldn't be" : "Point isn't in 3D simplex (triangle), but should be") << "\n";
            print("Point", points[i]);
        }
    }

    coords = {{1., 4.}, {6., 3.}, {7.1, 6.5}};
    points = {{1., 4.}, {2., 4.2}, {3.5, 3.5}, {0.2, 0.2}, {7.1-1e-12, 6.5}};
    in_simplex = {true, true, true, false, false};
    for (size_t i = 0; i < points.size(); i++)
    {
        if (is_in_simplex(points[i], coords) != in_simplex[i])
        {
            error_output++;
            std::cout << (!in_simplex[i] ? "Point is in 3D simplex (triangle), but shouldn't be" : "Point isn't in 3D simplex (triangle), but should be") << "\n";
            print("Point", points[i]);
        }
    }

    // Is in tetrahedron? (4D simplex/tetrahedron)
    coords = {{0., 0., 0.}, {0., 0., 1.}, {0., 1., 0.}, {1., 0., 0.}};
    points = {{0.2, 0.2, 0.2}, {0.25, 0.25, 0.25}, {0., 1., 0.}, {1., 1., 0.}, {0., 1., 1.}, {0.5, 0.5, 0.5}, {1., 1., 1.}};
    in_simplex = {true, true, true, false, false, false, false};
    for (size_t i = 0; i < points.size(); i++)
    {
        if (is_in_simplex(points[i], coords) != in_simplex[i])
        {
            error_output++;
            std::cout << (!in_simplex[i] ? "Point is in 4D simplex (tetrahedron), but shouldn't be" : "Point isn't in 4D simplex (tetrahedron), but should be") << "\n";
            print("Point", points[i]);
        }
    }

    if (error_output > 0)
	{
		print("Errors occurred in test_simplex()", error_output);
	}
	else
	{
		print("No errors occurred in test_simplex()", error_output);
	}
    return error_output;
}

int test_cubic_roots()
{
    int error_output = 0;

    double a2{ -0.92455459054338163 }, a1{ 0.27900309322530026 }, a0{ -0.027600384437325549 };
    a2 = -1;
    a1 = 0.1356257109;
    a0 = -0.006523118931;
    a2 = -0.9222039261;
    a1 = 0.2834866938;
    a0 = -0.02904806022;
    std::vector<std::complex<double>> Z_analytical = cubic_roots_analytical(a2, a1, a0);
    print("Z", Z_analytical);
    std::vector<std::complex<double>> Z_iterative = cubic_roots_iterative(a2, a1, a0, 1e-14);
    print("Z", Z_iterative);

    return error_output;
}

double fn_pointer(std::function<double(double, double)> fun, double a, double b) { return fun(a, b); }
double add(double a, double b) { return a+b; }
int test_function_pass()
{
    int error_output = 0;
    auto f = std::bind(&add, std::placeholders::_1, std::placeholders::_2);
    error_output += (fn_pointer(f, 1., 1.5) == 2.5) ? 0 : 1;
    error_output += (fn_pointer(f, -1., 3.2) == 2.2) ? 0 : 1;
    if (error_output > 0)
	{
		print("Errors occurred in test_function_pass()", error_output);
	}
	else
	{
		print("No errors occurred in test_function_pass()", error_output);
	}
    return error_output;
}

struct Polynomial
{
    std::vector<double> a, b;
    double c;

    Polynomial(std::vector<double> a_) : a(a_), b(a_), c(0.) { }
    Polynomial(std::vector<double> a_, std::vector<double> b_, double c_) : a(a_), b(b_), c(c_) { }

    double evaluate(double x)
    {
        std::vector<double> d = (x < c) ? a : b;

        double y = 0.;
        int order = static_cast<int>(d.size()-1);
        for (int i = 0; i <= order; i++)
        {
            y += d[i] * std::pow(x, order-i);
        }
        return y;
    }
    double gradient(double x)
    {
        std::vector<double> d = (x < c) ? a : b;

        double dy = 0.;
        int order = static_cast<int>(d.size()-1);
        for (int i = 0; i <= order-1; i++)
        {
            dy += (order-i) * d[i] * std::pow(x, order-i-1);
        }
        return dy;
    }
};
int test_rootfinding()
{
    int error_output = 0;

    std::vector<std::pair<Polynomial, double>> references = {
        // Linear polynomial
        {Polynomial({3., 1.5}), -0.5},

        // Quadratic polynomial
        {Polynomial({2., -3., -1.5}), 1.895643924},

        // Cubic polynomial
        {Polynomial({1., 3., 2., 6.}), -3.},
        {Polynomial({1., 1., -5., 3.}), -3.},

        // Quartic polynomial
        {Polynomial({1., 0., 3., 0., -4}), 1.},

        // Quintic polynomial
        {Polynomial({1., 0., 2., 0., -0.5, 3.}), -1.044499064},

        // Discontinuous polynomials
        {Polynomial({3., 1.5}, {3., 1.5}, 0.), -0.5},
        {Polynomial({3., -0.5}, {3., 1.5}, 0.), 0.},
        {Polynomial({1., -2.5}, {1., -3., 3.}, 0.), 0.},
        {Polynomial({1., -2.5}, {1., -3., 3.}, 2.), 2.},
    };

    for (auto ref: references)
    {
        auto f = std::bind(&Polynomial::evaluate, ref.first, std::placeholders::_1);
        double root{ 0. }, a{ -5. }, b{ 5. };
        RootFinding rootfinding;
        int output = rootfinding.brent(f, root, a, b, 1e-12, 1e-14);
        root = rootfinding.getx();

        if (std::fabs(root - ref.second) > 1e-8 || output > 0)
        {
            print("Bisection root finding not successful", std::vector<double>{root, ref.second});
            error_output++;
        }
    }

    for (auto ref: references)
    {
        auto f = std::bind(&Polynomial::evaluate, ref.first, std::placeholders::_1);
        auto df = std::bind(&Polynomial::gradient, ref.first, std::placeholders::_1);
        double root{ 0. }, a{ -5. }, b{ 5. };
        RootFinding rootfinding;
        int output = rootfinding.bisection_newton(f, df, root, a, b, 1e-12, 1e-14);
        root = rootfinding.getx();

        if (std::fabs(root - ref.second) > 1e-8 || output > 0)
        {
            print("Bisection-Newton root finding not successful", std::vector<double>{root, ref.second});
            error_output++;
        }
    }

    for (auto ref: references)
    {
        auto f = std::bind(&Polynomial::evaluate, ref.first, std::placeholders::_1);
        double root{ 0. }, a{ -5. }, b{ 5. };
        RootFinding rootfinding;
        int output = rootfinding.brent(f, root, a, b, 1e-12, 1e-14);
        root = rootfinding.getx();

        if (std::fabs(root - ref.second) > 1e-8 || output > 0)
        {
            print("Brent root finding not successful", std::vector<double>{root, ref.second});
            error_output++;
        }
    }

    for (auto ref: references)
    {
        auto f = std::bind(&Polynomial::evaluate, ref.first, std::placeholders::_1);
        auto df = std::bind(&Polynomial::gradient, ref.first, std::placeholders::_1);
        double root{ 0. }, a{ -5. }, b{ 5. };
        RootFinding rootfinding;
        int output = rootfinding.brent_newton(f, df, root, a, b, 1e-12, 1e-14);
        root = rootfinding.getx();

        if (std::fabs(root - ref.second) > 1e-8 || output > 0)
        {
            print("Brent-Newton root finding not successful", std::vector<double>{root, ref.second});
            error_output++;
        }
    }

    if (error_output > 0)
	{
		print("Errors occurred in test_rootfinding()", error_output);
	}
	else
	{
		print("No errors occurred in test_rootfinding()", error_output);
	}
    return error_output;
}

int test_line_search()
{
    int error_output = 0;

    // Test the line search method considering a minimization of f(x) problem
    // Ha = -b | H = grad(b) ; b = grad(f) ; a = delta(x)

    double lamb = 1.;
    bool check;

    // Application to polynomial functions
    // f(x,y) = x1^4 + x1^2 + x2^2 | x = [ x1 x2 ]

    LineSearch  linesearch{2};
    // Add gradient b
    Eigen::VectorXd b = Eigen::VectorXd::Zero(2);
    b << -6., -2.;
    // Add step a
    Eigen::VectorXd a = Eigen::VectorXd::Zero(2);
    a << -3., -1.;

    // New values of independent variable vector x_old before calculation
    std::vector<double> x = {-2., 0.};
    // Old values of independent variable vector x_old before calculation
    std::vector<double> x_old = {1., 1.};
    // Objective function f(x), new value
    double f = 20.;
    // Objective function f(x), old value
    double f_old = 3.;

    check = linesearch.init(lamb, x_old, f_old, b, a, 1.0);

    while(f > f_old && check)
    {
        if(f > f_old)
        {
            check = linesearch.process(x, f);
            lamb = linesearch.get_alam();

            for(size_t i = 0; i < x.size(); i++)
            {
                x[i] = x_old[i] + lamb*a[i];
            }
            f = std::pow(x[0],4) + std::pow(x[0],2) + std::pow(x[1],2);
        }
        else
        {
            break;
        }
    }

    // Reduced the objective function ?
    check = (f_old - f > 0.);

    if (check == false)
	{
        error_output += 1;
		print("Errors occurred in test_line_search()", error_output);
	}
	else
	{
		print("No errors occurred in test_line_search()", error_output);
	}
    return error_output;
}