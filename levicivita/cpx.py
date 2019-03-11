"""A class for calculating with infinite and infinitesimal complex numbers"""

import cmath
import numbers

# testing only
import unittest

from . import *

_UNARY_NAMES = '''exp log log10 sqrt 
                  sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh 
                  phase polar
                  isnan isinf isfinite'''.split()

__all__ = ('LeviCivitaComplex', 'epsilon', 'eps', 'd', 'ε', 'L',
           'log', 'isclose', 'st', 'change_terms',
           'e', 'nan', 'nanj', 'pi', 'tau', 'inf', 'infj') + tuple(_UNARY_NAMES)

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

for name in _UNARY_NAMES:
    exec(f"""
def {name}(x, **kw):
    try:
        return x.{name}(**kw)
    except AttributeError:
        return cmath.{name}(x, **kw)
""", globals())

del name

# cmath.log(0) raises, but cmath.log(e) returns -inf+nanj
def log(x, base=None):
    b = cmath.e if base is None else base
    try:
        return x.log(b)
    except AttributeError:
        pass
    if isinstance(base, LeviCivitaBase):
        return type(base)(x).log(base)
    if base is None:
        return cmath.log(x)
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

def rect(r, phi):
    return LeviCivitaComplex.rect(r, phi)

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

from .real import LeviCivitaFloat as _REAL
LeviCivitaComplex._REAL = _REAL
del _REAL
