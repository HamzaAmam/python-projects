"""Transparent scalar root solvers using only Python's standard library."""

from dataclasses import dataclass
import math
from numbers import Real
from typing import Callable

Function = Callable[[float], float]


class ConvergenceError(RuntimeError):
    """The method stopped without satisfying its convergence rule."""


@dataclass(frozen=True)
class Iteration:
    """One evaluated approximation; iteration 0 means an initial point."""

    number: int
    x: float
    fx: float


@dataclass(frozen=True)
class RootResult:
    root: float
    residual: float  # Absolute value of f(root), in the equation's output units.
    iterations: int
    reason: str
    history: tuple[Iteration, ...]


def _finite(value: Real, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number.")
    try:
        number = float(value)
    except OverflowError as exc:
        raise ValueError(f"{name} is too large for floating-point arithmetic.") from exc
    if not math.isfinite(number):
        raise ValueError(f"{name} must be finite.")
    return number


def _controls(x_tol: float, max_iter: int) -> float:
    tolerance = _finite(x_tol, "x_tol")
    if tolerance <= 0:
        raise ValueError("x_tol must be positive.")
    if isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter < 1:
        raise ValueError("max_iter must be a positive integer.")
    return tolerance


def _evaluate(function: Function, x: float, name: str = "f") -> float:
    try:
        return _finite(function(x), f"{name}({x})")
    except (ArithmeticError, ValueError, TypeError) as exc:
        raise ValueError(f"Cannot evaluate {name} at x={x}: {exc}") from exc


def _result(history: list[Iteration], reason: str) -> RootResult:
    last = history[-1]
    return RootResult(last.x, abs(last.fx), last.number, reason, tuple(history))


def bisection(
    f: Function, a: float, b: float, *, x_tol: float = 1e-8, max_iter: int = 100
) -> RootResult:
    """Halve a sign-changing bracket until its half-width is <= x_tol.

    Caller must supply a continuous function on [a, b]. A sign change alone
    cannot distinguish a root from a discontinuity. x_tol has the units of x.
    """
    if not callable(f):
        raise ValueError("f must be callable.")
    x_tol = _controls(x_tol, max_iter)
    a, b = _finite(a, "a"), _finite(b, "b")
    if a >= b:
        raise ValueError("Require a < b.")
    fa, fb = _evaluate(f, a), _evaluate(f, b)
    for x, fx in ((a, fa), (b, fb)):
        if fx == 0:
            return _result([Iteration(0, x, fx)], "exact endpoint root")
    # Compare signs directly; multiplying fa * fb can overflow or underflow.
    if (fa > 0) == (fb > 0):
        raise ValueError("Endpoints must have opposite function signs.")

    history = []
    for number in range(1, max_iter + 1):
        # Half each endpoint first to avoid overflowing b - a.
        midpoint = a / 2 + b / 2
        fm = _evaluate(f, midpoint)
        history.append(Iteration(number, midpoint, fm))
        half_width = b / 2 - a / 2
        if fm == 0:
            return _result(history, "exact function zero")
        if half_width <= x_tol:
            return _result(history, "bracket half-width <= x_tol")
        if midpoint == a or midpoint == b:
            raise ConvergenceError("Bisection reached floating-point resolution; increase x_tol.")
        if (fa > 0) != (fm > 0):
            b = midpoint
        else:
            a, fa = midpoint, fm
    raise ConvergenceError(f"Bisection did not converge in {max_iter} iterations.")


def newton_raphson(
    f: Function,
    derivative: Function,
    x0: float,
    *,
    x_tol: float = 1e-8,
    f_tol: float = 1e-8,
    max_iter: int = 100,
    derivative_tol: float = 1e-14,
) -> RootResult:
    """Apply x_next = x - f(x)/f'(x), recording every approximation.

    After the initial check, require both a small step and small residual
    (or an exact function zero). This is not a rigorous root-error bound.
    derivative_tol is an absolute threshold with units of f'(x).
    """
    if not callable(f) or not callable(derivative):
        raise ValueError("f and derivative must be callable.")
    x_tol = _controls(x_tol, max_iter)
    f_tol = _finite(f_tol, "f_tol")
    derivative_tol = _finite(derivative_tol, "derivative_tol")
    if f_tol <= 0 or derivative_tol <= 0:
        raise ValueError("f_tol and derivative_tol must be positive.")
    x = _finite(x0, "x0")
    fx = _evaluate(f, x)
    history = [Iteration(0, x, fx)]
    if abs(fx) <= f_tol:
        return _result(history, "initial residual <= f_tol")

    for number in range(1, max_iter + 1):
        slope = _evaluate(derivative, x, "derivative")
        if abs(slope) <= derivative_tol:
            raise ConvergenceError("Derivative is too small; choose another guess or use bisection.")
        next_x = x - fx / slope
        if not math.isfinite(next_x):
            raise ConvergenceError("Newton step became non-finite; choose another guess.")
        next_fx = _evaluate(f, next_x)
        history.append(Iteration(number, next_x, next_fx))
        if next_fx == 0:
            return _result(history, "exact function zero")
        if abs(next_x - x) <= x_tol and abs(next_fx) <= f_tol:
            return _result(history, "step <= x_tol and residual <= f_tol")
        if next_x == x:
            raise ConvergenceError("Newton stagnated with residual above f_tol.")
        x, fx = next_x, next_fx
    raise ConvergenceError(f"Newton-Raphson did not converge in {max_iter} iterations.")
