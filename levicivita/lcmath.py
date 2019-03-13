"""Module replacing mathematical functions on complex values (from cmath)
with functions on Levi-Civita complex (including infinite and infinitesimal).
"""

import cmath

from . import *

# NOTE: There is no cmath.gamma. We could implement lcmath.gamma anyway,
# just without the fallthrough?

_UNARY_NAMES = '''exp log10 sqrt
                  sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh 
                  phase polar
                  isnan isinf isfinite'''.split()

__all__ = ('LeviCivitaComplex', 'epsilon', 'eps', 'd', 'ε', 'L',
           'log', 'isclose', 'st', 'derivative', 'change_terms',
           'e', 'nan', 'nanj', 'pi', 'tau', 'inf', 'infj') + tuple(_UNARY_NAMES)

for name in _UNARY_NAMES:
    exec(f"""
def {name}(x, **kw):
    try:
        return x.{name}(**kw)
    except AttributeError:
        pass
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
        pass
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

# TODO: should this really override parent?
L = LeviCivitaComplex
