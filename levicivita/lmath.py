"""Module replacing mathematical functions on real values (from math)
with functions on Levi-Civita real (including infinite and infinitesimal).
"""

import math

from . import *

# math functions not covered (as of Python 3.7):
#   factorial, gamma, lgamma
#   frexp, ldexp
#   fsum
#   gcd
#   remainder
#   expm1, log1p
#   erf, erfc

_UNARY_NAMES = '''exp log10 log2 sqrt hypot fabs modf
                  degrees radians
                  sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh 
                  isnan isinf isfinite'''.split()

_BINARY_NAMES = '''atan2 copysign hypot'''.split()

_DUNDER_NAMES = '''ceil floor round trunc'''.split()

_EXTRA_NAMES = tuple(_UNARY_NAMES + _BINARY_NAMES + _DUNDER_NAMES)

__all__ = ('LeviCivitaFloat', 'epsilon', 'eps', 'd', 'ε', 'L',
           'log', 'isclose', 'st', 'derivative', 'change_terms',
           'e', 'nan', 'pi', 'tau', 'inf') + _EXTRA_NAMES

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

# TODO: should this really override parent?
L = LeviCivitaFloat
