# Python Engineering Projects

## Project 01: Nonlinear Equation Solver

An educational Python project connecting numerical methods with mechanical engineering. It implements **bisection** and **Newton-Raphson** from scratch, with iteration histories, input validation, automated tests, and two worked engineering calculations.

**Status:** Implemented and locally verified with 25 passing automated tests. Uses only Python's standard library. Examples are illustrative models with checkable analytical answers; they are not experimental validation or a certified engineering design tool.

### Installation and quick start

Install Python 3.10 or newer and Git, then run these commands in a terminal:

```bash
git clone https://github.com/HamzaAmam/python-projects.git
cd python-projects
python --version
python -m examples.engineering
python -m examples.engineering --trace
python -m unittest discover -s tests -v
```

There are **no third-party dependencies and no pip installation step**. Run the commands from the repository root so Python can find the local packages. On Windows, substitute `py` for `python` if needed. On macOS/Linux, `python3` may be the command instead.

Without Git: use GitHub's **Code → Download ZIP**, extract it, open a terminal inside the extracted `python-projects-main` folder, and run the same Python commands.

### Engineering examples

| Model | Equation solved | Inputs | Expected answer |
| --- | --- | --- | --- |
| Nonlinear spring | `1000*x + 10000*x**3 - 110 = 0` | Linear stiffness 1000 N/m, cubic stiffness 10000 N/m³, load 110 N | Deflection `x = 0.100000000 m` |
| Spherical tank | `pi*h**2*(3*R-h)/3 - V = 0` | Radius `R = 1 m`, volume `V = 5*pi/24 m³` | Liquid depth `h = 0.500000000 m` |

Each example runs both solvers and prints the root, absolute residual, iteration count, and stopping reason. `--trace` prints every evaluated approximation. Last digits and iteration counts can vary with floating-point arithmetic. The spherical example happens to hit its root exactly in bisection; that does not mean bisection is generally faster.

### Solve your own equation

For `x² - 2 = 0`, save the following as `my_calculation.py` in the repository root and run `python my_calculation.py`:

```python
from root_finding import ConvergenceError, bisection, newton_raphson


def f(x):
    return x**2 - 2


def df(x):
    return 2*x


try:
    bracketed = bisection(f, 1, 2, x_tol=1e-10)
    tangent = newton_raphson(f, df, 1, x_tol=1e-10, f_tol=1e-10)
    print("Bisection:", bracketed.root)
    print("Newton:", tangent.root)
    for step in tangent.history:
        print(step.number, step.x, step.fx)
except (ValueError, ConvergenceError) as error:
    print("Calculation failed:", error)
```

Both answers should be close to `1.41421356237`. Supply Python functions directly; there is no expression parser or `eval()`.

### Method selection and stopping rules

| Method | Required inputs | Successful stopping rule | Main limitation |
| --- | --- | --- | --- |
| Bisection | Continuous `f`, ordered endpoints with opposite signs (or an exact endpoint root) | Exact function zero, or current bracket half-width ≤ `x_tol` | Cannot check continuity; may miss roots without a sign change |
| Newton-Raphson | `f`, analytical derivative, initial guess | Initial residual ≤ `f_tol`; subsequently exact zero or both step ≤ `x_tol` and residual ≤ `f_tol` | Can diverge, cycle, stagnate, or leave the physical domain |

`x_tol` is an **absolute** tolerance in the units of the unknown. `f_tol` is an absolute tolerance in the equation's output units. Defaults are `1e-8`. For the spring, those units are meters and newtons respectively. A small residual is not automatically a small error in x, especially near a flat function. Bisection uses a bracket-based x-error criterion and does not promise a residual below `f_tol`.

Newton's `derivative_tol` defaults to `1e-14` and rejects absolute slopes at or below that threshold. Its units are those of the derivative; rescaling an equation may require adjusting this threshold and `f_tol`. Both methods default to `max_iter=100`.

Successful calls return a frozen `RootResult` with `root`, absolute `residual`, `iterations`, `reason`, and a tuple of `Iteration(number, x, fx)` records. Newton history includes iteration 0 (the initial guess). Bisection history records midpoints, except an endpoint root is recorded at iteration 0. `iterations` counts updates/midpoints, not total function evaluations. Failures raise an exception rather than returning a misleading successful result.

### Error handling

- `ValueError`: invalid numbers, nonpositive tolerances, invalid iteration limits, noncallable inputs, a missing sign-changing bracket, or a function/derivative that returns an invalid value or raises a numerical/domain error.
- `ConvergenceError`: iteration limit reached, Newton derivative too small, nonfinite Newton step, or floating-point stagnation.

A solver cannot verify that your equation, derivative, units, or physical assumptions are correct. Newton does not impose bounds or fall back to bisection. For the tank, always check `0 ≤ h ≤ 2R` and avoid Newton guesses at the zero-slope endpoints. The example uses an interior initial guess and its physical range is tested.

### Project layout

```text
root_finding/
    __init__.py         Public solver imports
    solvers.py          Algorithms, validation, result records
examples/
    __init__.py
    engineering.py      Spring and tank models; command-line trace
tests/
    test_solvers.py     25 standard-library unittest tests
docs/
    learning_guide.md    Mathematics and a code walkthrough
```

### Tests and learning path

The tests check known roots, engineering answers, derivatives against finite differences, histories, endpoint handling, invalid inputs, domain errors, nonfinite arithmetic, iteration limits, a Newton cycle, and stagnation. Tests run locally using `unittest`; no hosted CI workflow is configured.

Read [the learning guide](docs/learning_guide.md), run the trace, then work through its exercises. Start with `examples/engineering.py` to understand the equations; read `root_finding/solvers.py` next to see the loops. Suggested extensions in the guide are **not implemented**.
