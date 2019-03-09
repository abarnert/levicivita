"""A class for calculating with infinite and infinitesimal real numbers"""

import math
import functools
import numbers

# testing only
import unittest

from . import *
from .cpx import LeviCivitaComplex

_UNARY_NAMES = '''exp log log10 sqrt 
                  sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh 
                  isnan isinf isfinite'''.split()

_DUNDER_NAMES = '''ceil floor round trunc'''.split()

__all__ = ('LeviCivitaFloat', 'epsilon', 'eps', 'd', 'ε', 'L',
           'log', 'isclose', 'st', 'change_terms',
           'e', 'nan', 'pi', 'tau', 'inf') + tuple(_UNARY_NAMES) + tuple(_DUNDER_NAMES)

# TODO: It would be nice if LeviCivitaFloat(2) * 1j gave you a
# LeviCivitaComplex(2j) instead of a TypeError.

@functools.total_ordering
class LeviCivitaFloat(LeviCivitaBase, numbers.Real):
    _MATH = math
    _TYPE = float
    _ABSTYPE = numbers.Real

    def __lt__(self, other):
        if not isinstance(other, numbers.Real):
            raise TypeError(f"'<' not supported between instances of '{type(self).__name__}' and '{type(other).__name__}'")
        if not isinstance(other, LeviCivitaBase):
            other = type(self)(other)

        # IEEE inf overwhelms everything else
        if self.isnan() or other.isnan():
            return False
        if math.isinf(self.front) or math.isinf(other.front):
            return self.front < other.front

        # handle zero and negative numbers
        if self.front == 0 or other.front == 0:
            return self.front < other.front
        if self.front < 0 and other.front > 0:
            return True
        if self.front > 0 and other.front < 0:
            return False
        if self.front < 0 and other.front < 0:
            return -other < -self

        # TODO: Maybe just convert from front/leading/series to series and
        # compare lexicographically?
        
        if self.leading < other.leading:
            return False
        if other.leading < self.leading:
            return True

        if self.front < other.front:
            return True
        if other.front < self.front:
            return False

        return self - other < 0

    # ABC checks before total_ordering can help...
    def __le__(self, other):
        return self == other or self < other
    
    def __ceil__(self):
        if self.isinf():
            raise OverflowError('cannot convert infinity to integer')
        return math.ceil(self.st())

    def __floor__(self):
        if self.isinf():
            raise OverflowError('cannot convert infinity to integer')
        return math.floor(self.st())

    def __round__(self, ndigits=None):
        if self.isinf():
            raise OverflowError('cannot convert infinity to integer')
        return round(self.st(), ndigits)
    
    def __trunc__(self):
        if self.isinf():
            raise OverflowError('cannot convert infinity to integer')
        return math.trunc(self.st())
    
    def __float__(self):
        return float(self.st())

    def __floordiv__(self, other):
        raise NotImplementedError # TODO

    def __rfloordiv__(self, other):
        raise NotImplementedError # TODO

    def __mod__(self, other):
        raise NotImplementedError # TODO

    def __rmod__(self, other):
        raise NotImplementedError # TODO

    # TODO: There must be a way to refactor this nicely, but the need
    #       to realify intermediate values makes it tricky...
    
    @classmethod
    def _realify(cls, c, test=False):
        if test and not math.isclose(c.front.imag, 0, abs_tol=1e-15):
            raise ValueError('math domain error')
        return cls(c.front.real, c.leading,
                   tuple((q, a.real) for (q, a) in c.series))

    def exp(self):
        return self._realify(LeviCivitaComplex(self).exp(), test=True)
    
    def log(self, base=math.e):
        return self._realify(LeviCivitaComplex(self).log(base), test=True)

    def log10(self):
        return self.log(10)

for name in 'cos sin tan acos asin atan cosh sinh tanh acosh asinh atanh'.split():
    exec(f'''
def {name}(self):
    return self._realify(LeviCivitaComplex(self).{name}())
LeviCivitaFloat.{name} = {name}
del {name}
    ''', globals())
    
for name in _UNARY_NAMES:
    exec(f"""
def {name}(x, **kw):
    try:
        return x.{name}(**kw)
    except AttributeError:
        return math.{name}(x, **kw)
""", globals())

# TODO: optional argument on round!
for name in _DUNDER_NAMES:
    exec(f"""
def {name}(x, **kw):
    try:
        return x.__{name}__(**kw)
    except AttributeError:
        return math.{name}(x, **kw)
""", globals())
    
del name

def log(x, base=math.e):
    try:
        return x.log(base)
    except AttributeError:
        pass
    if isinstance(base, LeviCivitaBase):
        return type(base)(x).log(base)
    return math.log(x, base)

def isclose(x, y, *, rel_tol=1e-09, abs_tol=0.0):
    try:
        return x.isclose(y, rel_tol=rel_tol, abs_tol=abs_tol)
    except AttributeError:
        pass
    try:
        return y.isclose(x, rel_tol=rel_tol, abs_tol=abs_tol)
    except AttributeError:
        pass
    return math.isclose(x, y, rel_tol=rel_tol, abs_tol=abs_tol)

def st(x):
    try:
        return x.st()
    except AttributeError:
        return float(x)

epsilon = eps = d = ε = LeviCivitaFloat(1.0, 1)
# TODO: should these all be LeviCivitaFloat(n)?
e = math.e
nan = math.nan
pi = math.pi
tau = math.tau
# TODO: should these be 1/ε instead of inf? or should there be another constant?
inf = math.inf

L = LeviCivitaFloat # TODO: move to test/__main__ only?
