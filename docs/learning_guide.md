# Understanding the solver

## 1. Turn an engineering question into f(x) = 0

A root is an input that makes a function zero. Move every term of your equation to one side so the solver can measure the imbalance.

For the spring, assume a static load and the illustrative force law:

```text
F = k*x + alpha*x^3
f(x) = k*x + alpha*x^3 - F
f'(x) = k + 3*alpha*x^2
```

Here x is deflection in meters, k is 1000 N/m, alpha is 10000 N/m³, and F is 110 N. This model assumes the given positive cubic stiffness describes the spring; it does not model dynamics, damping, yielding, or measured material behavior.

At x = 0.1 m, the linear term is 100 N and the cubic term is 10 N. Their sum is 110 N, so f(0.1) = 0 (up to computer rounding). Both stiffness terms are increasing for positive x, giving a unique positive equilibrium.

For the tank, the liquid forms a spherical cap measured from the bottom of a sphere:

```text
V(h) = pi*h^2*(3R - h)/3
f(h) = V(h) - V_target
f'(h) = pi*h*(2R - h)
```

Use R = 1 m and target volume 5*pi/24 m³ (about 0.654498 m³). Substitute h = 0.5 m to verify this exact target volume. A physical spherical tank has 0 <= h <= 2R. The formula assumes an ideal sphere and ignores internal hardware. Other algebraic roots outside this depth range are not physically valid.

## 2. Bisection: keep cutting a valid bracket in half

If f is continuous and the endpoints have opposite signs, the intermediate value theorem guarantees at least one root between them. Evaluate the midpoint, then keep the half whose endpoints still have opposite signs.

For f(x) = x² - 2 on [1, 2]:

| Iteration | a | b | Midpoint | f(midpoint) | Next bracket |
| --- | --- | --- | --- | --- | --- |
| 1 | 1 | 2 | 1.5 | 0.25 | [1, 1.5] |
| 2 | 1 | 1.5 | 1.25 | -0.4375 | [1.25, 1.5] |
| 3 | 1.25 | 1.5 | 1.375 | -0.109375 | [1.375, 1.5] |

At the nth evaluated midpoint, the maximum distance to a root in the retained bracket is `(b_initial-a_initial)/2**n`, in exact arithmetic. For [1, 2] and x_tol = 1e-8, 27 midpoint evaluations suffice unless an exact zero is found sooner. Floating-point limits still apply.

The code compares signs instead of multiplying function values, avoiding unnecessary overflow and underflow. It computes the midpoint as `a/2 + b/2` to avoid overflowing `b-a` for very large opposite-sign endpoints.

A trap: f(x) = 1/x changes sign across zero but is discontinuous there. There is no root. The caller must check continuity; a sign change is not proof by itself. An even-multiplicity root, such as f(x) = (x-1)², may have no sign change and cannot generally be bracketed this way.

## 3. Newton-Raphson: follow a tangent to the axis

Near a guess x, approximate f with a straight line:

```text
f(x + delta) approximately equals f(x) + f'(x)*delta
0 = f(x) + f'(x)*delta
delta = -f(x)/f'(x)
x_next = x - f(x)/f'(x)
```

For x² - 2, start at x = 1. The next guesses are 1.5, 1.416666667, and 1.414215686. Newton often converges quadratically close to a simple root when the function is sufficiently smooth and the guess is suitable. That behavior is conditional, not a promise.

A zero or very small derivative makes the tangent step undefined or enormous. A poor initial guess can cycle: x³ - 2x + 2 starting at 0 jumps to 1 and back to 0 indefinitely. One automated test verifies that this reaches the iteration limit instead of reporting success.

## 4. Error, residual, and stopping are different things

- Root error is the distance from an approximation to the true root. Usually the true root is unknown.
- Residual is `abs(f(x))`. It measures imbalance in the equation's output units.
- Step size is `abs(x_next - x)`. A tiny step can reflect stagnation, not accuracy.

Bisection stops on a sufficiently narrow bracket or an exact function zero. Under the continuity assumption, the bracket gives an x-error bound, subject to rounding.

Newton accepts an initial point if its residual is already small. After an update, it requires both a small step and a small residual, unless f evaluates exactly to zero. These checks are practical convergence criteria, not a rigorous bound on root error. A nearly flat function can have a small residual far from its root. An exact floating-point zero also depends on how your function was evaluated.

All tolerances are absolute. Scale them for the units and accuracy your model warrants. Requesting 30 decimal places from ordinary Python floats is not a clever accuracy setting; it is a request the arithmetic cannot fulfill.

## 5. Read the Python in this order

1. In `examples/engineering.py`, each `def` names a function that maps one input to one output. Constants carry unit comments.
2. In `solvers.py`, the `Function` type hint describes a callable; it does not execute it. The solver receives a function rather than hardcoding one engineering equation.
3. `_finite` and `_controls` reject invalid values before the loop. `_evaluate` adds the failing x value to numerical error messages.
4. A `for` loop caps the work at `max_iter`. `if` statements implement stopping and bracket selection. `raise` reports a failed calculation to the caller.
5. `@dataclass(frozen=True)` creates readable result records that cannot be casually modified. `history` becomes a tuple when returned.
6. The `if __name__ == '__main__'` guard runs the example only when invoked, not when tests import its functions.
7. In the tests, assertions compare results with independently known answers. Central differences provide a separate numerical check on the analytical derivatives.

## 6. Exercises (not implemented)

1. Trace x² - 2 by hand for three steps, then compare with `history`.
2. Change the spring load to 200 N and solve again. Check force balance by substituting your computed x back into the original equation. Update the expected answer and tests when changing constants.
3. Run the Newton cycle above and catch `ConvergenceError`. Explain why increasing max_iter will not fix that cycle.
4. Try a different tank volume between zero and 4*pi*R³/3. Choose a physical bracket and verify the returned depth and volume.
5. Add a secant method or a plot of residual versus iteration as a future extension, with its own tests.

Before presenting this as portfolio work, be able to explain the sign check, derive the Newton update, state each tolerance's units, and demonstrate a failed calculation. Running code is the first step; understanding why it works is the point.
