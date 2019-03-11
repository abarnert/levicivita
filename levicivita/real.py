"""A class for calculating with infinite and infinitesimal real numbers"""

import math
import functools
import numbers

from . import *

# math functions not covered (as of 3.7):
#   factorial, gamma, lgamma
#   frexp, ldexp
#   fsum
#   gcd
#   remainder
#   expm1, log1p
#   erf, erfc

_UNARY_NAMES = '''exp log log10 log2 sqrt hypot fabs modf
                  degrees radians
                  sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh 
                  isnan isinf isfinite'''.split()

_BINARY_NAMES = '''atan2 copysign hypot'''.split()

_DUNDER_NAMES = '''ceil floor round trunc'''.split()

_EXTRA_NAMES = tuple(_UNARY_NAMES + _BINARY_NAMES + _DUNDER_NAMES)

__all__ = ('LeviCivitaFloat', 'epsilon', 'eps', 'd', 'ε', 'L',
           'log', 'isclose', 'st', 'change_terms',
           'e', 'nan', 'pi', 'tau', 'inf') + _EXTRA_NAMES

# TODO: It would be nice if LeviCivitaFloat(2) * 1j gave you a
# _COMPLEX(2j) instead of a TypeError.

@functools.total_ordering
class LeviCivitaFloat(LeviCivitaBase, numbers.Real):
    _MATH = math
    _TYPE = float
    _ABSTYPE = numbers.Real
    _COMPLEX = None

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

    def sign(self):
        if self < 0:
            return -1
        elif self > 0:
            return 1
        else:
            # handle -0.0, but only if there's no infinitesimal
            return math.copysign(1, self.front)

    @staticmethod
    def _sign(value):
        # -0.0+ε is positive, so only fall back to math.copysign for
        # actual 0
        if value < 0:
            return -1
        elif value > 0:
            return 1
        else:
            if isinstance(value, LeviCivitaFloat):
                return math.copysign(1, value.front)
            else:
                return math.copysign(1, value)
        
    def copysign(self, other):
        ssign, osign = self._sign(self), self._sign(other)
        flip = ssign*osign
        front = self.front*flip if self.front else math.copysign(0.0, osign)
        return type(self)(front, self.leading, self.series)
    
    def fabs(self):
        return abs(self)

    def hypot(self, other):
        return (self*self + other*other).sqrt()

    def degrees(self):
        return self * math.degrees(1)

    def radians(self):
        return self * math.radians(1)
    
    def modf(self):
        intpart = trunc(self)
        return self-intpart, intpart
    
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
        return floor(self / other)

    def __rfloordiv__(self, other):
        return floor(other / self)

    def __mod__(self, other):
        return self - self//other

    def __rmod__(self, other):
        return other - floor(other/self)

    def conjugate(self):
        return self
    
    @property
    def real(self):
        return self

    @property
    def imag(self):
        raise type(self)()

    # TODO: There must be a way to refactor this nicely, but the need
    #       to realify intermediate values makes it tricky... Although
    #       even realifying intermediates doesn't keep the errors in
    #       check, so maybe it isn't worth doing?
    
    @classmethod
    def _realify(cls, c, test=False):
        if test and not math.isclose(c.front.imag, 0, abs_tol=1e-15):
            raise ValueError('math domain error')
        return c.real

    def log(self, base=math.e):
        if self <= 0:
            raise ValueError('math domain error')
        return super().log(base)

    def log10(self):
        return self.log(10)

    def log2(self):
        return self.log(2)

    def atan2(self, other):
        if other > 0:
            return (self/other).atan()
        elif other < 0:
            if self >= 0:
                return (self/other).atan() + math.pi
            else:
                return (self/other).atan() - math.pi
        else:
            if self >= 0:
                return type(self)(math.pi/2)
            elif self < 0:
                return type(self)(-math.pi/2)
            else:
                # math.atan2(0, 0) is pi/2, not NaN
                return type(self)()

for name in _UNARY_NAMES:
    exec(f"""
def {name}(x, **kw):
    try:
        return x.{name}(**kw)
    except AttributeError:
        pass
    return math.{name}(x, **kw)
""", globals())

# TODO: optional argument on round!
for name in _DUNDER_NAMES:
    exec(f"""
def {name}(x, **kw):
    try:
        return x.__{name}__(**kw)
    except AttributeError:
        pass
    return math.{name}(x, **kw)
""", globals())

for name in _BINARY_NAMES:
    exec(f"""
def {name}(x, y, **kw):
    try:
        return x.{name}(y, **kw)
    except AttributeError:
        pass
    if isinstance(y, LeviCivitaFloat):
        return type(y)(x).{name}(y, **kw)
    if isinstance(y, LeviCivitaBase):
        raise TypeError(f"can't convert {{type(y).__name__}} to LeviCivitaFloat")
    return math.{name}(x, y, **kw)
""", globals())

del name

def log(x, y=math.e):
    try:
        return x.log(y)
    except AttributeError:
        pass
    if isinstance(y, LeviCivitaFloat):
        return type(y)(x).log(y)
    if isinstance(y, LeviCivitaBase):
        raise TypeError(f"can't convert {{type(y).__name__}} to LeviCivitaFloat")
    return math.log(x, y)
    
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
        pass
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

from .cpx import LeviCivitaComplex as _COMPLEX
LeviCivitaFloat._COMPLEX = _COMPLEX
del _COMPLEX
