"""Two illustrative engineering models, with SI units throughout."""

import argparse
import math

from root_finding import RootResult, bisection, newton_raphson


def spring_force_balance(x: float) -> float:
    """k*x + alpha*x**3 - F = 0; x in m, residual in N."""
    k = 1000.0  # N/m
    alpha = 10000.0  # N/m^3; positive cubic stiffness
    force = 110.0  # N
    return k * x + alpha * x**3 - force


def spring_derivative(x: float) -> float:
    """Derivative of force balance, in N/m."""
    return 1000.0 + 30000.0 * x**2


def tank_volume_balance(h: float) -> float:
    """Spherical-cap volume minus target; h in m, residual in m^3.

    Radius R=1 m. Physical depth is restricted to 0 <= h <= 2R.
    Target 5*pi/24 m^3 is chosen to give the checkable depth h=0.5 m.
    """
    radius = 1.0
    target_volume = 5.0 * math.pi / 24.0
    return math.pi * h**2 * (3.0 * radius - h) / 3.0 - target_volume


def tank_derivative(h: float) -> float:
    """dV/dh = pi*h*(2R-h), in m^2, with R=1 m."""
    return math.pi * h * (2.0 - h)


def print_result(label: str, result: RootResult, residual_unit: str, trace: bool) -> None:
    print(f"{label}: {result.root:.9f} m | residual {result.residual:.3e} {residual_unit}"
          f" | {result.iterations} iterations")
    print(f"  Stop: {result.reason}")
    if trace:
        print("  iteration         x (m)                 f(x)")
        for step in result.history:
            print(f"  {step.number:9d} {step.x:20.12g} {step.fx:20.12g}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trace", action="store_true", help="show each approximation")
    args = parser.parse_args()
    print("Nonlinear spring: 1000*x + 10000*x^3 = 110 N; expected x = 0.1 m")
    print_result("Bisection", bisection(spring_force_balance, 0.0, 0.3), "N", args.trace)
    print_result("Newton", newton_raphson(spring_force_balance, spring_derivative, 0.05), "N", args.trace)
    print("\nSpherical tank: R = 1 m, V = 5*pi/24 m^3; expected h = 0.5 m")
    print_result("Bisection", bisection(tank_volume_balance, 0.0, 2.0), "m^3", args.trace)
    print_result("Newton", newton_raphson(tank_volume_balance, tank_derivative, 0.8), "m^3", args.trace)


if __name__ == "__main__":
    main()
