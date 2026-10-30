"""Run from repository root: python -m unittest discover -s tests -v."""
import math
import unittest

from examples.engineering import (
    spring_derivative, spring_force_balance, tank_derivative, tank_volume_balance,
)
from root_finding import ConvergenceError, bisection, newton_raphson


class BisectionTests(unittest.TestCase):
    def test_square_root_and_history(self):
        result = bisection(lambda x: x*x - 2, 1, 2, x_tol=1e-10)
        self.assertLessEqual(abs(result.root - math.sqrt(2)), 1e-10)
        self.assertEqual(result.iterations, len(result.history))
        self.assertEqual(result.root, result.history[-1].x)
        self.assertEqual(result.residual, abs(result.root**2 - 2))

    def test_endpoint_roots(self):
        for a, b in ((2, 3), (1, 2)):
            with self.subTest(a=a, b=b):
                result = bisection(lambda x: x - 2, a, b)
                self.assertEqual(result.root, 2)
                self.assertEqual(result.iterations, 0)

    def test_midpoint_root(self):
        self.assertEqual(bisection(lambda x: x, -1, 1).root, 0)

    def test_negative_root(self):
        self.assertAlmostEqual(bisection(lambda x: x*x - 2, -2, -1).root, -math.sqrt(2), places=7)

    def test_bad_intervals(self):
        for a, b in ((2, 1), (1, 1), (-1, 1)):
            with self.subTest(a=a, b=b):
                with self.assertRaises(ValueError):
                    bisection(lambda x: x*x + 1, a, b)

    def test_iteration_limit(self):
        with self.assertRaises(ConvergenceError):
            bisection(lambda x: x*x - 2, 1, 2, max_iter=1)

    def test_large_bracket_midpoint(self):
        self.assertEqual(bisection(lambda x: x, -1e308, 1e308).root, 0)

    def test_tiny_function_signs_do_not_underflow(self):
        with self.assertRaises(ValueError):
            bisection(lambda x: 1e-300, 0, 1)

    def test_floating_point_stagnation(self):
        a, b = 1.0, math.nextafter(1.0, 2.0)
        with self.assertRaises(ConvergenceError):
            bisection(lambda x: -1 if x == a else 1, a, b, x_tol=1e-30)


class NewtonTests(unittest.TestCase):
    def test_square_root_and_history(self):
        result = newton_raphson(lambda x: x*x - 2, lambda x: 2*x, 1)
        self.assertAlmostEqual(result.root, math.sqrt(2), places=10)
        self.assertLessEqual(result.residual, 1e-8)
        self.assertEqual(len(result.history), result.iterations + 1)
        self.assertEqual(result.history[0].x, 1)

    def test_initial_root_needs_no_derivative(self):
        self.assertEqual(newton_raphson(lambda x: x*x, lambda x: 0, 0).iterations, 0)

    def test_zero_derivative(self):
        with self.assertRaises(ConvergenceError):
            newton_raphson(lambda x: x*x + 1, lambda x: 2*x, 0)

    def test_cycle_hits_iteration_limit(self):
        # x=0 -> 1 -> 0: Newton is not guaranteed to converge.
        with self.assertRaises(ConvergenceError):
            newton_raphson(lambda x: x**3 - 2*x + 2, lambda x: 3*x*x - 2, 0, max_iter=10)

    def test_small_step_with_large_residual_is_not_success(self):
        with self.assertRaises(ConvergenceError):
            newton_raphson(lambda x: 1, lambda x: 1e30, 1)

    def test_nonfinite_step(self):
        with self.assertRaises(ConvergenceError):
            newton_raphson(lambda x: 1e308, lambda x: 1e-10, 1)

    def test_iteration_limit(self):
        with self.assertRaises(ConvergenceError):
            newton_raphson(lambda x: x*x - 2, lambda x: 2*x, 1, max_iter=1)


class ValidationTests(unittest.TestCase):
    def test_invalid_common_controls(self):
        for options in ({"x_tol": 0}, {"x_tol": -1}, {"x_tol": math.nan},
                        {"x_tol": True}, {"max_iter": 0}, {"max_iter": 1.5},
                        {"max_iter": True}):
            for solve in (
                lambda: bisection(lambda x: x, -1, 1, **options),
                lambda: newton_raphson(lambda x: x, lambda x: 1, 1, **options),
            ):
                with self.subTest(options=options):
                    with self.assertRaises(ValueError):
                        solve()

    def test_invalid_numbers(self):
        for value in (math.nan, math.inf, -math.inf, True, "1", 1j):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    bisection(lambda x: x, value, 2)
                with self.assertRaises(ValueError):
                    newton_raphson(lambda x: x, lambda x: 1, value)

    def test_invalid_newton_controls(self):
        for options in ({"f_tol": 0}, {"f_tol": math.inf}, {"derivative_tol": -1},
                        {"derivative_tol": math.nan}):
            with self.subTest(options=options):
                with self.assertRaises(ValueError):
                    newton_raphson(lambda x: x, lambda x: 1, 1, **options)

    def test_nonfinite_function_and_derivative(self):
        for value in (math.nan, math.inf, 1j, True):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    bisection(lambda x: value, -1, 1)
                with self.assertRaises(ValueError):
                    newton_raphson(lambda x: x, lambda x: value, 1)

    def test_domain_error(self):
        with self.assertRaisesRegex(ValueError, "Cannot evaluate"):
            bisection(math.log, -1, 2)

    def test_noncallable(self):
        with self.assertRaises(ValueError):
            bisection(4, 0, 1)
        with self.assertRaises(ValueError):
            newton_raphson(lambda x: x, None, 1)


class EngineeringTests(unittest.TestCase):
    def test_spring_matches_known_deflection(self):
        # Independent check: 1000*0.1 + 10000*0.1^3 = 110 N.
        for result in (bisection(spring_force_balance, 0, 0.3),
                       newton_raphson(spring_force_balance, spring_derivative, 0.05)):
            self.assertAlmostEqual(result.root, 0.1, places=7)
            self.assertLess(abs(spring_force_balance(result.root)), 2e-5)

    def test_tank_matches_known_depth_in_physical_domain(self):
        # Independent check: pi*(1/2)^2*(3-1/2)/3 = 5*pi/24.
        for result in (bisection(tank_volume_balance, 0, 2),
                       newton_raphson(tank_volume_balance, tank_derivative, 0.8)):
            self.assertAlmostEqual(result.root, 0.5, places=7)
            self.assertTrue(0 <= result.root <= 2)

    def test_analytical_derivatives_against_central_difference(self):
        for f, derivative, x in ((spring_force_balance, spring_derivative, 0.1),
                                 (tank_volume_balance, tank_derivative, 0.5)):
            with self.subTest(x=x):
                delta = 1e-6
                numerical = (f(x + delta) - f(x - delta)) / (2*delta)
                self.assertAlmostEqual(numerical, derivative(x), delta=1e-6)


if __name__ == "__main__":
    unittest.main()
