"""Classes for calculating with infinite and infinitesimal real numbers"""

import abc
import cmath
import collections
import fractions
import functools
import itertools
import math
import numbers
import operator
import sys
import types

# debug only
import inspect

__version__ = '0.0.9'

__all__ = ('LeviCivitaBase', 'LeviCivitaFloat', 'LeviCivitaComplex',
           'L', 'change_terms',
           'epsilon', 'eps', 'd', 'ε',
           'derivative')

def _debugmethod(func):
    return func

    @functools.wraps(func)
    def meth(self, *args):
        depth = len(inspect.stack())
        print(f'{depth*">"} {self!r}.{meth.__name__}{args}')
        result = func(self, *args)
        print(f'{depth*"<"} {self!r}.{meth.__name__}{args} = {result!r}')
        return result
    return meth

class LeviCivitaBase(abc.ABC):
    # form is front*d^leading*(sum a_q*d^q), where all q's are >0
    __slots__ = ('front', 'leading', 'series')

    def __new__(cls, front=None, leading=None, series=None):
        """Construct a Levi-Civita number"""
        self = super(LeviCivitaBase, cls).__new__(cls)

        # Zero
        if front is None:
            self.front, self.leading = cls._TYPE(0), 0
            self.series = ((0, 1.0),)
            return self

        # Copy
        if leading is None:
            # TODO: construct from string? Decimal?
            try:
                self.front, self.leading = cls._TYPE(front.front), front.leading
                self.series = front.series
                if (not isinstance(front, cls._ABSTYPE) or
                    not all(isinstance(a, cls._ABSTYPE)
                            for (q, a) in front.series)):
                    raise TypeError(f"can't convert {type(front).__name__} to {cls.__name__}")
                return self
            except AttributeError:
                if isinstance(front, cls._ABSTYPE):
                    self.front, self.leading = cls._TYPE(front), 0
                    self.series = ((0, 1.0),)
                    return self
                raise TypeError(f"can't convert {type(front).__name__} to {cls.__name__}")

        if not isinstance(front, cls._ABSTYPE):
            raise TypeError(f'front must be {cls._ABSTYPE.__name__}, not {front}')
        # TODO: handle math.inf, -math.inf, math.nan, cmath.nan?
        if not isinstance(leading, numbers.Rational):
            try:
                leading = fractions.Fraction(leading)
            except ValueError:
                raise TypeError(f'front must be rational, not {leading}')
        if not series:
            series = ((0, 1.0),)
        # TODO: don't need to sort here; this is just type-checking...
        try:
            series = tuple(sorted((tuple(term) for term in series),
                                  key=operator.itemgetter(0)))
        except TypeError:
            raise TypeError(f'series must be an iterable of pairs, not {series}')
        for term in series:
            if not isinstance(term[0], numbers.Rational):
                # TODO: change to Fraction?
                raise TypeError(f'first element of each term must be rational, not {series[0]}')
            if not isinstance(term[1], numbers.Complex):
                raise TypeError(f'second element of each term must be real or complex, not {series[1]}')
        front, leading, series = cls._normalize(front, leading, series)
        self.front, self.leading, self.series = front, leading, series
        return self

    @staticmethod
    def _collect(series):
        last_q = last_aq = None
        for q, aq in series:
            if q == last_q:
                last_aq += aq
            else:
                if last_aq:
                    yield (last_q, last_aq)
                last_q, last_aq = q, aq
        if last_aq:
            yield (last_q, last_aq)

    @classmethod
    @_debugmethod
    def _normalize(cls, front, leading, series):
        #print(f'normalizing {front}, {leading}, {series}')
        # TODO: Do we really need to sort, instead of just building a Counter?
        series = tuple(cls._collect(sorted(series, key=operator.itemgetter(0))))
        if not series:
            return 0.0, 0, ((0, 1),)
        series = series[:cls._TERMS]
        #print(f'  collected {series}')
        k = series[0][1]
        front *= k
        series = tuple((q, aq/k) for (q, aq) in series)
        #print(f'  mul-normalized to {front}, {leading}, {series}')
        p = series[0][0]
        if p:
            leading += p
            series = tuple((q-p, aq) for (q, aq) in series)
            #print(f'  add-normalized to {front}, {leading}, {series}')
        return cls._TYPE(front), leading, series

    def __repr__(self):
        return f'{type(self).__name__}({self.front}, {self.leading}, {self.series})'

    def __str__(self):
        if self.isnan():
            return 'nan'
        if not self.front:
            return '0'
        parts = []
        for i, (q, a) in enumerate(self.series[:self._DISPLAY_TERMS]):
            power = q + self.leading
            coeff = a * self.front
            if coeff:
                if coeff == -1 and power:
                    coeff = '-'
                elif coeff == 1 and power:
                    coeff = '+' if i else ''
                else:
                    coeff = f'{coeff:+g}' if i else f'{coeff:g}'
                    if '/' in coeff or '+' in coeff[1:] or '-' in coeff[1:]:
                        if coeff[0] in '+-':
                            coeff = coeff[0] + '(' + coeff[1:] + ')'
                        else:
                            coeff = '(' + coeff + ')'
                    if power:
                        coeff += '*'
                if power == 0:
                    part = coeff
                elif power == 1:
                    part = f'{coeff}ε'
                else:
                    part = f'{coeff}ε**{power}'
                parts.append(part)
        return ''.join(parts)

    def st(self):
        if self.leading < 0:
            raise OverflowError('cannot take standard part of infinite')
        elif self.leading > 0:
            return 0
        else:
            return self.front

    def __bool__(self):
        return bool(self.front)

    # TODO: verify NaN handling in all arithmetic (also inf?)

    # Hook to allow coerce, e.g., `LeviCivitaReal(1) + 1j` to
    # LeviCivitaComplex (because __radd__ won't help there).
    def _coerce_binop(meth):
        @functools.wraps(meth)
        def wrapper(self, other, *args, **kw):
            if (not isinstance(other, LeviCivitaBase)
                and isinstance(self, numbers.Real)
                and isinstance(other, numbers.Complex)
                and not isinstance(other, numbers.Real)):
                try:
                    self = self._COMPLEX(self)
                except AttributeError:
                    pass
            return meth(self, other, *args, **kw)
        return wrapper

    def __abs__(self):
        # TODO: we probably need to abs the coefficients in the series?
        return type(self)(abs(self.front), self.leading, self.series)

    @_debugmethod
    def __neg__(self):
        return type(self)(-self.front, self.leading, self.series)

    @_debugmethod
    def __pos__(self):
        return self

    @_debugmethod
    @_coerce_binop
    def __add__(self, other):
        if not isinstance(other, type(self)):
            try:
                other = type(self)(other)
            except TypeError:
                return NotImplemented
        if not self.front:
            return other # not just optimization; avoids divide by zero
        h = other.leading - self.leading
        ff = other.front / self.front
        series = tuple((term[0]+h, term[1]*ff) for term in other.series)
        return type(self)(self.front, self.leading, self.series + series)

    __radd__ = __add__

    @_debugmethod
    @_coerce_binop
    def __sub__(self, other):
        return self + -other

    @_debugmethod
    @_coerce_binop
    def __rsub__(self, other):
        return -self + other

    @_coerce_binop
    def __eq__(self, other):
        if not isinstance(other, type(self)):
            try:
                other = type(self)(other)
            except TypeError:
                return NotImplemented
        # TODO: inf incorrectly returns L(0, 1) != L(0, 0) which implies
        #       that maybe the right behavior is more complicated than this?
        if not self.front and not other.front:
            return True
        return ((self.front, self.leading, self.series) ==
                (other.front, other.leading, other.series))

    @_coerce_binop
    def __ne__(self, other):
        return not self == other

    @_coerce_binop
    def __mul__(self, other):
        if not self:
            return self
        if not isinstance(other, type(self)):
            try:
                other = type(self)(other)
            except TypeError:
                return NotImplemented
        if not other:
            return other
        front = self.front * other.front
        leading = self.leading + other.leading
        series = tuple((x[0]+y[0], x[1]*y[1])
                       for x in self.series
                       for y in other.series)
        return type(self)(front, leading, series)

    __rmul__ = __mul__

    def eps_part(self):
        return type(self)(1, 0, self.series) - 1

    # expand a Taylor series
    def expand(self, coefficients):
        s = type(self)()
        pow = type(self)(1.0) # = self**0
        for coefficient in itertools.islice(coefficients, self._TERMS):
            s += pow * coefficient
            pow *= self
        return s

    # expand a Taylor series with inverse coefficients
    # (except that 0 is 0, not 1/0)
    def expand_inv(self, inv_coefficients):
        s = type(self)()
        pow = type(self)(1.0) # = self**0
        for inv_coefficient in itertools.islice(inv_coefficients, self._TERMS):
            if inv_coefficient:
                s += pow / inv_coefficient
            pow *= self
        return s

    def inv(self):
        # reduce it to inverting 1/(1-e)
        z = type(self)(1/self.front, -self.leading)
        return z * (-self.eps_part()).expand(itertools.repeat(1))

    @_coerce_binop
    def __truediv__(self, other):
        if not isinstance(other, type(self)):
            try:
                other = type(self)(other)
            except TypeError:
                return NotImplemented
        return self * other.inv()

    @_coerce_binop
    def __rtruediv__(self, other):
        return self.inv() * other

    def sqrt(self):
        return self**fractions.Fraction(1, 2)

    def _intpow(self, p):
        if p == 0:
            return type(self)(1.0, 0)
        if p == 1:
            return self
        if p == 2:
            return self*self
        if not self.front:
            return type(self)() if p > 0 else type(self)(self._MATH.nan)
        if p < 0:
            return self.inv() ** -p
        r = self ** (p//2)
        r = r * r
        if p % 2:
            r = r * self
        return r

    @_coerce_binop
    def __pow__(self, p):
        if isinstance(p, numbers.Integral):
            return self._intpow(p)
        if isinstance(p, LeviCivitaBase):
            return (self.log() * p).exp()
        if not self.leading:
            return type(self)(self._MATH.exp(self.log() * p))
        if not isinstance(p, numbers.Rational):
            try:
                p = fractions.Fraction(p)
            except ValueError:
                raise TypeError(f'{type(self).__name__} can only take rational powers')
        # TODO: Real front often gives a tiny imaginary part. Is there
        # a way to either avoid this, or correct for it after the fact?
        front = self.front ** p
        leading = self.leading / p.denominator
        z = type(self)(front, leading)
        def f(i, u):
            return u*(p-1+i)/1 if i else 1
        taylor = self._generate_taylor(f)
        return z * self.eps_part().expand(taylor)

    def __rpow__(self, other):
        try:
            other = type(self)(other)
        except TypeError:
            return NotImplemented
        return other ** self

    def conjugate(self):
        raise NotImplementedError # TODO

    @property
    def real(self):
        raise NotImplementedError # TODO

    @property
    def imag(self):
        raise NotImplementedError # TODO

    def __complex__(self):
        return complex(self.st())

    @classmethod
    def _isnan(cls, value):
        if isinstance(value, LeviCivitaBase):
            return value.isnan()
        return cls._MATH.isnan(value)

    def isnan(self):
        return self._isnan(self.front)

    @classmethod
    def _isinf(cls, value):
        if isinstance(value, LeviCivitaBase):
            return value.isinf()
        return cls._MATH.isinf(value)

    def isinf(self):
        # TODO: Should this really be asking "is IEEE infinity or an
        #       infinite L-C", or only the former (in which case we'd
        #       want a separate method/function for the latter)?
        return self._isinf(self.front) or self.leading < 0

    def isfinite(self):
        return not self.isinf()

    def isclose(self, other, *, rel_tol=1e-09, abs_tol=0.0):
        if self._isnan(self) or self._isnan(other):
            return False
        if self._isinf(self) or self._isinf(other):
            other = type(self)(other)
            return (self.leading == other.leading and
                    cmath.isclose(self.front, other.front, rel_tol=rel_tol))
        if isinstance(other, LeviCivitaBase):
            other = other.st()
        # TODO: abs?
        return cmath.isclose(self.st(), other, rel_tol=rel_tol, abs_tol=abs_tol)

    def exp(self):
        # This diverges for infinite values, to infinities too large
        # to approximate in a finite Levi-Civita series. Which makes
        # sense once you think about it.
        if self.leading < 0:
            raise OverflowError('cannot exponentiate infinite')
        # TODO: this seems to be the one place where 5/10 terms is
        # nowhere near enough...
        return self.expand_inv(math.factorial(i) for i in itertools.count())

    def log(self, base=math.e):
        # TODO: is this correct?
        if not self.front:
            raise ValueError('cannot take log of 0')
        if self.leading < 0:
            raise ValueError('cannot take log of infinite')
        elif self.leading > 0:
            raise ValueError('cannot take log of infinitesimal')
        if isinstance(base, LeviCivitaBase):
            return self.log() / base.log()
        if base != math.e:
            return self.log() / cmath.log(base)
        series = (i * -1**(i-1) for i in itertools.count())
        return self._MATH.log(self.front) + self.eps_part().expand(series)

    def log2(self):
        return self.log(2)

    def log10(self):
        return self.log(10)

    # Some of the inverse trig/hyp identities use coefficients based
    # on double factorials (1*3*5 / 2*4*6, etc.). Writing a simple
    # `factorial2` function and then writing each of the identities in
    # closed form is probably the cleanest way to do this, but since
    # they all happen to use exactly the same sequence (modulo signs),
    # I'm just generating that sequence here.
    @staticmethod
    def _doublefactorialseries():
        coeff = fractions.Fraction(1)
        for i in itertools.count():
            if not i%2:
                yield 0
            else:
                yield coeff / i
                coeff *= fractions.Fraction(i, i+1)

    def cos(self):
        # 1, 0, -2, 0, 24, 0, -720, 0, ...
        series = (math.factorial(i) * (1,0,-1,0)[i%4]
                  for i in itertools.count())
        return self.expand_inv(series)

    def sin(self):
        # 0, 1, 0, -6, 0, 120, 0, -5040, ...
        series = (math.factorial(i) * (0,1,0,-1)[i%4]
                  for i in itertools.count())
        return self.expand_inv(series)

    def tan(self):
        # TODO: expand directly?
        return self.sin() / self.cos()

    def acos(self):
        return self._MATH.pi/2 - self.asin()

    def asin(self):
        # 0, 1, 0, 1/6, 0, 3/40, 0, 5/112, 0, ...
        series = (n * (0,1,0,1)[i%4]
                  for i, n in enumerate(self._doublefactorialseries()))
        return self.expand(series)

    def atan(self):
        # 0, 1, 0, -3, 0, 5, 0, -7, ...
        series = (i * (0,1,0,-1)[i%4] for i in itertools.count())
        return self.expand_inv(series)

    def cosh(self):
        # 1, 0, 2, 0, 24, 0, 720, 0, ...
        series = (math.factorial(i) * (1,0,1,0)[i%4]
                  for i in itertools.count())
        return self.expand_inv(series)

    def sinh(self):
        # 0, 1, 0, 6, 0, 120, 0, 5040, ...
        series = (math.factorial(i) * (0,1,0,1)[i%4]
                  for i in itertools.count())
        return self.expand_inv(series)

    def tanh(self):
        # TODO: expand directly?
        return self.sinh() / self.cosh()

    def acosh(self):
        return (self + (self**2 - 1).sqrt()).log()

    def asinh(self):
        # 0, 1, 0, -1/6, 0, 3/40, 0, -5/112, 0, ...
        #series = (n * (0,1,0,-1)[i%4] for n in self._doublefactorialseries())
        #return self.expand(series)
        return (self + (self**2 + 1).sqrt()).log()

    def atanh(self):
        return ((1 + self) / (1 - self)).log() / 2

    @classmethod
    def _generate_taylor(cls, f, *, terms=None):
        if terms is None:
            terms = cls._TERMS
        t = []
        # TODO: reduce/accumulate?
        value = None
        for i in range(cls._TERMS):
            value = f(i, value)
            t.append(value)
        return t

    # TODO: Context? Class properties?
    @classmethod
    def change_terms(cls, terms=5):
        cls._DISPLAY_TERMS = terms
        cls._TERMS = 2*terms
        cls._TAYLOR_EXP = cls._generate_taylor(lambda i, l: l/i if i else 1)
        cls._TAYLOR_LOG = cls._generate_taylor(lambda i, l: -(l/abs(l))/i if i>1 else i)

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

class LeviCivitaComplex(LeviCivitaBase, numbers.Complex):
    _MATH = cmath
    _TYPE = complex
    _ABSTYPE = numbers.Complex
    _REAL = None

    def __abs__(self):
        return super().__abs__().real

    def conjugate(self):
        return type(self)(self.front.conjugate(), self.leading,
                          ((q, a.conjugate()) for (q, a) in self.series))

    @property
    def real(self):
        return self._REAL(self.front.real, self.leading,
                          ((q, a.real) for (q, a) in self.series))

    @property
    def imag(self):
        return self._REAL(self.front.imag, self.leading,
                          ((q, a.imag) for (q, a) in self.series))

    def phase(self):
        return self.imag.atan2(self.real)

    @classmethod
    def rect(cls, r, phi):
        x = cls(r * cos(phi))
        y = cls(r * sin(phi))
        return x + y*j

    def polar(self):
        return self.abs(), self.phase()

def L(*args, **kw):
    def _complex(val):
        return (isinstance(arg, numbers.Complex) and
                not isinstance(arg, numbers.Real))
    args = list(args)
    for i, arg in enumerate(args):
        if _complex(arg):
            return LeviCivitaComplex(*args, **kw)
        if isinstance(arg, collections.Iterable):
            args[i] = arg = tuple(arg)
            if any(_complex(part) for part in arg):
                return LeviCivitaComplex(*args, **kw)
    return LeviCivitaFloat(*args, **kw)

LeviCivitaFloat._COMPLEX = LeviCivitaComplex
LeviCivitaComplex._REAL = LeviCivitaFloat

change_terms = LeviCivitaBase.change_terms
change_terms()

epsilon = eps = d = ε = LeviCivitaFloat(1.0, 1)

from .diff import derivative
