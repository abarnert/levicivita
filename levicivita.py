# A class for calculating with infinity and infinitesimals
# Initial version is based on the JavaScript "inf" calculator
#   https://github.com/bcrowell/inf/blob/master/LeviCivita.js
#   http://www.lightandmatter.com/calc/inf/
# Also see
#   https://en.wikipedia.org/wiki/Levi-Civita_field
#   http://www2.physics.umanitoba.ca/u/khodr/Publications/RS-Overview-offprints.pdf
#
# The Levi-Civita field is an extension of the real numbers with one
# added infinitesimal element eps, and then other infinitesimal and
# infinite numbers are constructed from there. Each member is a
# series: sum_q(a_q*(eps**q)), where each a_q is a real number and
# each q is a rational number. The support (the set of indices where
# a_q != 0) is a left-finite set: for any q, there are only finitely
# many terms less than it with nonzero a_q. Ordering is lexicographical
# ordering of the components. A real number is a number where the only
# nonzero element (if any) is a_0.
#
# Levi-Civita infinitesimals are not as powerful as hyperreals; they
# can only be used for calculus on functions that have Taylor series
# expansions. However, those Taylor series can be carried out in IEEE
# floats, keeping track of only the first few terms (by default, the
# inf calculator keeps 8 infinitesimal terms and displays 4), which
# allows computations with infinitesimal numbers that are as accurate
# as IEEE floating point to be automated. This is good enough for
# automatic differentiation, which is the main use for it.
#
# In the case of underflow of the noninfinitesimal part, you have to
# chose whether to raise, drop to 0, drop to eps, drop to a very large
# constant infinitesimal, drop to an arbitrary infinitesimal that's
# still larger than the infinitesimal part, etc. (And likewise for
# overflow.) I believe the inf calculator underflows to 0--or, in some
# cases, "1/inf", which is smaller than 1/eps. For example:
#   (10+eps)^-400 = (1/inf) + (-1/inf)eps + (1/inf)eps^2 + ...

#import sympy
#ε = eps = sympy.Symbol('ε')

import cmath
import fractions
import functools
import inspect
import math
import numbers
import operator
import sys
import types
import unittest

__version__ = '0.0.1'

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

# TODO: real as well as complex
#   includes adding comparison operators
class LeviCivitaComplex(numbers.Complex):
    # form is front*d^leading*(sum a_q*d^q), where all q's are >0
    __slots__ = ('front', 'leading', 'series')
    
    def __new__(cls, front=None, leading=None, series=None):
        """Construct a Levi-Civita number"""
        self = super(LeviCivitaComplex, cls).__new__(cls)

        # Zero
        if front is None:
            self.front, self.leading = 0.0, 0
            self.series = ((0, 1),)
            return self
        
        # Copy
        if leading is None:
            # TODO: construct from string? Decimal?
            # TODO: constructing from Fraction gives rational front,
            #       which is probably a bad idea
            try:
                self.front, self.leading = front.front, front.leading
                self.series = front.series
                return self
            except AttributeError:
                if isinstance(front, numbers.Complex):
                    self.front, self.leading = front, 0
                    self.series = ((0, 1),)
                    return self
                raise TypeError(f'single-argument {cls.__name__}() requires a LeviCivitaComplex, not {type(front)}')

        if not isinstance(front, numbers.Complex):
            raise TypeError(f'front must be real or complex, not {front}')
        # TODO: handle math.inf, -math.inf, math.nan, cmath.nan?
        if not isinstance(leading, numbers.Rational):
            try:
                leading = fractions.Fraction(leading)
            except ValueError:
                raise TypeError(f'front must be rational, not {leading}')
        if not series:
            series = ((0, 1),)
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
        return front, leading, series

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

    # TODO: return LeviCivitaReal once we have such a thing
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
    def __sub__(self, other):
        return self + -other

    @_debugmethod    
    def __rsub__(self, other):
        return -self + other

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

    def __ne__(self, other):
        return not self == other

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
        return 1.0 - type(self)(1, 0, self.series)

    # expand a Taylor series
    def expand(self, coefficients):
        s = type(self)()
        pow = type(self)(1.0) # = self**0
        for term in coefficients:
            s += term * pow
            pow *= self
        return s

    # TODO: This is only useful for debugging, and not even that
    #       anymore now that repr works. It _might_ be useful to
    #       iterate the actual terms, rather than the internal
    #       series, but if so, that shouldn't be __iter__.
    #def __iter__(self):
    #    for term in self.series[:self._TERMS]:
    #        yield (term[0] + self.leading, term[1] * self.front)

    def inv(self):
        # reduce it to inverting 1/(1-e)
        z = type(self)(1/self.front, -self.leading)
        return z * (-self.eps_part()).expand(self._TAYLOR_INV)

    def __truediv__(self, other):
        if not isinstance(other, type(self)):
            try:
                other = type(self)(other)
            except TypeError:
                return NotImplemented
        return self * other.inv()

    def __rtruediv__(self, other):
        return self.inv() * other

    def __pow__(self, p):
        if isinstance(p, numbers.Integral):
            if p == 0:
                return type(self)(1.0, 0)
            if p == 1:
                return self
            if p == 2:
                return self*self
            if not self.front:
                return type(self)() if p > 0 else type(self)(math.nan)
            if p < 0:
                return self.inv() ** -p
            r = self ** (p//2)
            r = r * r
            if p % 2:
                r = r * self
            return r
        if not self.leading:
            return type(self)(cmath.exp(self.log() * p))
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
        # TODO: Maybe we should implement this in terms of exp and log
        # instead of a custom Taylor series?
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

    # TODO: Most of the functions below will give complex results even
    # with real input. We could force to real if input is real and output
    # has a tiny infinitesimal, or even force to real whenever input is
    # real. But it's probably better to create a separate Real class.
    
    def exp(self):
        # TODO: is this correct?
        if self.leading < 0:
            raise OverflowError('cannot exponentiate infinite')
        magf = abs(self.front)
        if magf == 1:
            return self.expand(self._TAYLOR_EXP)
        else:
            unitf = cmath.rect(1, cmath.phase(self.front))
            z = type(self)(unitf, self.leading, self.series)
            return z.expand(z._TAYLOR_EXP) ** magf
        
    def log(self, base=math.e):
        # TODO: is this correct?
        if not self.front:
            raise ValueError('cannot take log of 0')
        if self.leading < 0:
            raise ValueError('cannot take log of infinite')
        elif self.leading > 0:
            raise ValueError('cannot take log of infinitesimal')
        if base != math.e:
            return self.log() / cmath.log(base)
        return cmath.log(self.front) + self.eps_part().expand(self._TAYLOR_LOG)

    def log10(self):
        return self.log(10)
    
    def cos(self):
        u = (self*1j).exp()
        return (u + u.inv()) / 2

    def sin(self):
        u = (self*1j).exp()
        return (u - u.inv()) / 2j

    def tan(self):
        return self.sin() / self.cos()

    def acos(self):
        return -j * (self + (self**2 - 1).sqrt()).log()
    
    def asin(self):
        return -j * (self*1j + (1 - self**2).sqrt()).log()

    def atan(self):
        u = self*1j
        return -.5j * ((1 - u)/(1 + u)).log()
    
    def cosh(self):
        u = self.exp()
        return (u + u.inv()) / 2

    def sinh(self):
        u = self.exp()
        return (u - u.inv()) / 2j

    def tanh(self):
        return self.sinh() / self.cosh()

    def acosh(self):
        return (self + (self**2 - 1).sqrt()).log()

    def asinh(self):
        return (self + (self**2 + 1).sqrt()).log()        

    def atanh(self):
        return ((1 + self) / (1 - self)).log() / 2

    def sqrt(self):
        return self**fractions.Fraction(1, 2)
    
    def isnan(self):
        return cmath.isnan(self.front)

    def isinf(self):
        # TODO: Should this really be asking "is IEEE infinity or an
        #       infinite L-C", or only the former (in which case we'd
        #       want a separate method/function for the latter)?
        return cmath.isinf(self.front) or self.leading < 0

    def isfinite(self):
        return not self.isinf()

    def isclose(self, other, *, rel_tol=1e-09, abs_tol=0.0):
        if isnan(self) or isnan(other):
            return False
        if isinf(self) or isinf(other):
            other = type(self)(other)
            return (self.leading == other.leading and
                    cmath.isclose(self.front, other.front, rel_tol=rel_tol))
        return cmath.isclose(st(self), st(other),
                             rel_tol=rel_tol, abs_tol=abs_tol)

    def phase(self):
        raise NotImplementedError # TODO

    @classmethod
    def rect(cls, r, phi):
        raise NotImplementedError # TODO

    @classmethod
    def polar(self):
        raise NotImplementedError # TODO
        
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
        cls._TAYLOR_INV = cls._generate_taylor(lambda i, l: 1)
        cls._TAYLOR_EXP = cls._generate_taylor(lambda i, l: l/i if i else 1)
        cls._TAYLOR_LOG = cls._generate_taylor(lambda i, l: -(l/abs(l))/i if i>1 else i)

change_terms = LeviCivitaComplex.change_terms
change_terms()

def _makeunary(name):
    exec(f"""
def {name}(x, **kw):
    try:
        return x.{name}(**kw)
    except AttributeError:
        return cmath.{name}(x, **kw)
""", globals())

for name in 'exp log10 sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh sqrt isnan isinf isfinite'.split():
    _makeunary(name)

del _makeunary

def log(x, base=math.e):
    try:
        return x.log(base)
    except AttributeError:
        pass
    if isinstance(base, LeviCivitaComplex):
        return type(base)(x).log(base)
    return cmath.log(x, base)

def isclose(x, y, *, rel_tol=1e-09, abs_tol=0.0):
    try:
        return x.isclose(y, rel_tol=rel_tol, abs_tol=abs_tol)
    except AttributeError:
        pass
    try:
        return y.isclose(x, rel_tol=rel_tol, abs_tol=abs_tol)
    except AttributeError:
        pass
    return cmath.isclose(x, y, rel_tol=rel_tol, abs_tol=abs_tol)

def st(x):
    try:
        return x.st()
    except AttributeError:
        return complex(x)

epsilon = eps = d = ε = LeviCivitaComplex(1.0, 1)
# TODO: should these all be LeviCivitaComplex(n)?
e = cmath.e
nan = cmath.nan
nanj = cmath.nanj
pi = cmath.pi
tau = cmath.tau
# TODO: should these be 1/ε instead of inf? or should there be another constant?
inf = cmath.inf
infj = cmath.infj

L = LeviCivitaComplex # TODO: move to test/__main__ only?

class TestLeviCivitaComplex(unittest.TestCase):
    def assertClose(self, first, second, msg=None, *,
                    rel_tol=1e-09, abs_tol=0.0):
        if isclose(first, second, rel_tol=rel_tol, abs_tol=abs_tol):
            return
        diff = abs(first - second)
        sr = unittest.util.safe_repr
        stdiff = diff if isinf(diff) else st(diff)
        standardMsg = f'{first} != {second} ({sr(first)} != {sr(second)}) within rel_tol {rel_tol}, abs_tol {abs_tol} ({stdiff} difference)'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)

    # TODO: This is a (bad) hack until we have reals
    def assertSmallerMagnitude(self, first, second, msg=None):
        mfirst, msecond = abs(L(first)), abs(L(second))
        if ((-mfirst.leading, mfirst.front, mfirst.series) <
            (-msecond.leading, msecond.front, msecond.series)):
            return
        sr = unittest.util.safe_repr
        standardMsg = f'|{first}| >= |{second}| (|{sr(mfirst)}| >= |{sr(msecond)}|)'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)
    
    def test_real(self):
        self.assertFalse(L(2)+L(-2))
        self.assertEqual(L(2)+L(-3), L(-1))
        self.assertEqual(L(2)+(-3), L(-1))
        self.assertEqual(L(2)-3, L(-1))
        self.assertEqual(2+L(-3), L(-1))
        self.assertEqual(2-L(3), L(-1))
        self.assertEqual(L(2)+L(-3), -1)
        with self.assertRaises(ZeroDivisionError):
            L(2)/L(0)

    def test_complex(self):
        self.assertEqual(L(1+2j) * L(2+1j), L(5j))
        self.assertClose(sqrt(L(-1)), 1j)
        self.assertClose(L(-1)**.5, 1j)
        self.assertClose(L(-1)**L(.5), 1j)
        with self.assertRaises(ZeroDivisionError):
            (2+3j)/L(0j)

    def test_eps(self):
        self.assertNotEqual(eps, 0)
        self.assertNotEqual(eps, L(0))
        self.assertTrue(eps)
        self.assertFalse(L(0))
        self.assertFalse(L(2)+eps-L(2)-eps)
        self.assertEqual(0, L(0))
        self.assertEqual(eps+eps, 2*eps)
        self.assertEqual(eps*2, 2*eps)
        self.assertEqual(L(2)*eps, 2*eps)
        self.assertEqual(eps*L(2), 2*eps)
        self.assertEqual(L(1)+eps, eps+L(1))
        self.assertEqual(L(1)+eps, L(1, 0, ((0, 1), (1, 1))))
        self.assertEqual((L(1)+eps).st(), 1)
        self.assertEqual(eps * (1/eps), 1)
        self.assertClose((1+2j+eps) * (2+1j+eps), 5j+(3+3j)*eps+eps**2)
        self.assertEqual(eps*eps, eps**2)
        self.assertEqual(eps*eps, pow(eps, 2))
        self.assertEqual(1/(1-eps), sum((-eps)**n for n in range(eps._TERMS)))
        self.assertEqual(exp(eps),
                         sum((1/math.factorial(i))*eps**i
                             for i in range(eps._TERMS)))
        self.assertClose(cos(sin(cos(sin(eps)))), .666367, rel_tol=1e-6)        
        self.assertEqual((1+eps)**pi, 1)
        # This is really testing internal implementaton details, so it's
        # probably not a great test...
        self.assertEqual(eps**pi, eps**(2**-48))
        self.assertClose(exp(log(2+3*d)), log(exp(2+3*d)), rel_tol=1e-6)
        with self.assertRaises(OverflowError):
            2**(1/eps)

    def test_compare(self):
        self.assertSmallerMagnitude(L(1), L(2))
        self.assertSmallerMagnitude(eps, 1e-300)
        self.assertSmallerMagnitude(eps, -1e-300)
        self.assertSmallerMagnitude(-eps, 1e-300)
        self.assertSmallerMagnitude(eps**2, eps)
        self.assertSmallerMagnitude(eps, eps*2)
        self.assertSmallerMagnitude(eps, sqrt(eps))
        
    def test_st(self):
        self.assertEqual(st(eps), 0)
        self.assertEqual(st(1+2j), 1+2j)
        self.assertEqual(st(L(1+2j)), 1+2j)
        self.assertEqual(st(1+2j+eps), 1+2j)
        self.assertEqual(st(1+2j+eps*2j), 1+2j)
        self.assertEqual(st((1+eps)*eps), 0)
        self.assertEqual(st((1+2j+eps) * (2+1j+eps)), 5j)
        
    def test_str(self):
        self.assertEqual(str(L(1)+eps), "1+ε")
        self.assertEqual(str(1-eps), "1-ε")
        self.assertEqual(str(1-eps**3), "1-ε**3")
        self.assertEqual(str((1-eps)**2), "1-2*ε+ε**2")
        self.assertEqual(str(1 + 1/eps), "ε**-1+1")
        self.assertEqual(str(1-eps**3-eps**-3), "-ε**-3+1-ε**3")
        
    # TODO add order comparisons once reals are added
    def test_underflow_overflow(self):
        # TODO: handle once we fix underflow
        min = sys.float_info.min
        self.assertNotEqual(min * eps, 0)
        #self.assertNotEqual(min * eps * min, 0)
        #self.assertNotEqual(min * (min * eps), 0)
        #self.assertNotEqual((min + eps)**2, 0)

        # TODO: handle once we fix overflow
        max_eps = sys.float_info.max + eps
        # currently gives inf+inf*eps
        # (inf calculator gives inf + nan/nan*eps + nan/nan*eps**2)
        #self.assert??((max+eps)**2, ???)
        
if __name__ == '__main__':
    unittest.main()
