"""Educational root-finding methods for engineering calculations."""

from .solvers import ConvergenceError, Iteration, RootResult, bisection, newton_raphson

__all__ = ["ConvergenceError", "Iteration", "RootResult", "bisection", "newton_raphson"]
