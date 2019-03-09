"""A class for calculating with infinite and infinitesimal complex numbers"""

import cmath
import numbers

# testing only
import unittest

from . import *

_UNARY_NAMES = '''exp log log10 sqrt 
                  sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh 
                  isnan isinf isfinite'''.split()

# TODO: phase, polar, rect

__all__ = ('LeviCivitaComplex', 'epsilon', 'eps', 'd', 'ε', 'L',
           'log', 'isclose', 'st', 'change_terms',
           'e', 'nan', 'nanj', 'pi', 'tau', 'inf', 'infj') + tuple(_UNARY_NAMES)

class LeviCivitaComplex(LeviCivitaBase, numbers.Complex):
    _MATH = cmath
    _TYPE = complex
    _ABSTYPE = numbers.Complex

    def phase(self):
        raise NotImplementedError # TODO

    @classmethod
    def rect(cls, r, phi):
        raise NotImplementedError # TODO

    @classmethod
    def polar(self):
        raise NotImplementedError # TODO

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
        if isinstance(base, LeviCivitaBase):
            return self.log() / base.log()
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

for name in 'exp log10 sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh sqrt isnan isinf isfinite'.split():
    exec(f"""
def {name}(x, **kw):
    try:
        return x.{name}(**kw)
    except AttributeError:
        return cmath.{name}(x, **kw)
""", globals())

del name

def log(x, base=math.e):
    try:
        return x.log(base)
    except AttributeError:
        pass
    if isinstance(base, LeviCivitaBase):
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

# TODO: phase, polar, rect
    
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
