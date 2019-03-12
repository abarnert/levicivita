import functools

from . import ε

__all__ = ('ε', 'derivative')

def derivative(func):
    """Return dfunc(x), the derivative of func(x).

    `func` must be defined in terms of operators, builtins like `abs`,
    and functions from `lmath` and `lcmath` only. If you use functions
    from `math` or `cmath`, or C extension modules like `numpy`, the
    derivative function will just raise a `TypeError` when called.
    """
    @functools.wraps(func)
    def dfunc(x):
        plus = (func(x+ε) - func(x+0*ε))/ε
        minus = (func(x-0*ε) - func(x-ε))/ε
        # TODO: isclose and st aren't quite right here. We want to
        # ignore rounding errors, even if they're on coefficients of
        # infinite terms. Notice that the tests (in test.py) on a simple
        # polynomial fail when given a huge or tiny complex number.
        # However, fixing this depends on deciding exactly what isclose
        # should mean.
        if not plus.isclose(minus):
            raise ValueError(f'{func.__name__} is not differentiable at {x}: {plus} != {minus}')
        return plus.st()
    dfunc.__name__ = 'd' + dfunc.__name__
    dfunc.__qualname__ = 'd' + dfunc.__qualname__
    return dfunc
