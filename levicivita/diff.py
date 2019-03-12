import numbers

from . import ε

__all__ = ('ε', 'derivative')

def derivative(func, x):
    """Return the value for the derivative of func at x.
 
    `LeviCivitaComplex` infinitesimals will be used if `lcmath` functions
    or a complex `x` are provided, `LeviCivitaFloat` otherwise.
    """
    plus = (func(x+ε) - func(x))/ε
    minus = (func(x) - func(x-ε))/ε
    if not plus.isclose(minus):
        raise ValueError(f'{func.__name__} is not differentiable at {x}: {plus} != {minus}')
    return plus.st()

if __name__ == '__main__':
    def test(func, dfunc, *xs):
        for x in xs:
            f = func(x)
            d1 = dfunc(x)
            try:
                d2 = derivative(func, x)
            except ValueError as e:
                d2 = str(e)
            print(f'  {x}: {d1} ?= {d2}')
        print()

    def func(x):
        return 2*x**3 + 3*x**2 + 4*x + 5
    def dfunc(x):
        return 6*x**2 + 6*x + 4
    print('2*x**3 + 3*x**2 + 4*x + 5')
    test(func, dfunc, 0, 0j, 3, 3+0j, 3j, 1e-100, 1e100)

    from levicivita.lmath import sin, cos, pi, L
    def func(x):
        return sin(x) * cos(x)
    def dfunc(x):
        return cos(x)**2 - sin(x)**2
    print('sin(x) * cos(x) (real)')
    test(func, dfunc, 0, L(1), L(pi/8), L(5*pi/32))

    from levicivita.lcmath import sin, cos, pi, L
    def func(x):
        return sin(x) * cos(x)
    def dfunc(x):
        return cos(x)**2 - sin(x)**2
    print('sin(x) * cos(x) (complex)')
    test(func, dfunc, 0j, L(1), L(1+0j), L(pi/8), L(5*pi/32), L(5j*pi/32))
