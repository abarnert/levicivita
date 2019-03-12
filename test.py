import cmath
import math
import pathlib
import sys
import unittest

from unittest.util import safe_repr as sr

from levicivita import lmath
from levicivita import lcmath

class _TestBaseLeviCivita(unittest.TestCase):
    def setUp(self):
        for name in 'ε pi exp log sin cos isclose isinf isnan sqrt st'.split():
            func = getattr(self._LCMOD, name)
            setattr(self, name, func)
        self.L = self._LCCLASS
        self.L.change_terms(10)
    
    def assertClose(self, first, second, msg=None, *,
                    rel_tol=1e-09, abs_tol=0.0):
        if self.isnan(first) and self.isnan(second):
            return
        if self.isclose(first, second, rel_tol=rel_tol, abs_tol=abs_tol):
            return
        diff = abs(first - second)
        stdiff = diff if self.isinf(diff) else self.st(diff)
        standardMsg = f'{first} != {second} ({sr(first)} != {sr(second)}) within rel_tol {rel_tol}, abs_tol {abs_tol} ({stdiff} difference)'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)
    
    def assertMathTest(self, funcname, value, *args, **kwargs):
        msg = kwargs.pop('msg', None)
        msg = (f'{msg}:{funcname}({value})' if msg else f'{funcname}({value})')
        # TODO: We do want to skip over functions we're not implementing,
        # like math.lgamma... but should there be some notice we're doing
        # so, to catch functions we didn't intend to skip?
        mathfunc = getattr(self._MATHMOD, funcname, None)
        lcfunc = getattr(self._LCMOD, funcname, None)
        if not mathfunc or not lcfunc:
            return
        try:
            expected = mathfunc(value)
        except Exception as e:
            ex = type(e)
        else:
            actual = lcfunc(value)
            # polar returns a pair of values
            if isinstance(expected, tuple):
                for e, a in zip(expected, actual):
                    self.assertClose(e, a, *args, msg=msg, **kwargs)
            else:
                self.assertClose(expected, actual, *args, msg=msg, **kwargs)
            return
        with self.assertRaises(ex, *args, msg=msg, **kwargs):
            lcfunc(value)

    def assertMathTestCase(self, line, *args, **kwargs):
        # We don't care about outval or flags, because we're testing
        # against math/cmath value or exception instead.
        testid, funcname, value, j_or_arrow, *rest = line.split()
        msg = kwargs.pop('msg', None)
        msg = f'{msg}:{testid}' if msg else testid
        value = self._MATHCLASS(value)
        if j_or_arrow != '->':
            value += float(j_or_arrow)*1j
        self.assertMathTest(funcname, value, *args, msg=msg, **kwargs)
        
    def assertMathTestCases(self, filename):
        # HACK because we're shadowing stdlib test package
        path = pathlib.Path(unittest.__file__).parent.parent / 'test'
        with open(path / filename) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('--'):
                    self.assertMathTestCase(line)

class TestLeviCivitaFloat(_TestBaseLeviCivita):
    _MATHMOD = math
    _MATHCLASS = float
    _LCMOD = lmath
    _LCCLASS = lmath.LeviCivitaFloat
    
    def assertLess(self, first, second, msg=None):
        if first < second:
            return
        sr = unittest.util.safe_repr
        standardMsg = f'{first} >= {second} ({sr(first)} >= {sr(second)})'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)
    
    def test_real(self):
        self.assertFalse(self.L(2)+self.L(-2))
        self.assertEqual(self.L(2)+self.L(-3), self.L(-1))
        self.assertEqual(self.L(2)+(-3), self.L(-1))
        self.assertEqual(self.L(2)-3, self.L(-1))
        self.assertEqual(2+self.L(-3), self.L(-1))
        self.assertEqual(2-self.L(3), self.L(-1))
        self.assertEqual(self.L(2)+self.L(-3), -1)
        with self.assertRaises(ZeroDivisionError):
            self.L(2)/self.L(0)

    def test_complex(self):
        with self.assertRaises(TypeError):
            self.L(lcmath.L(2))
        with self.assertRaises(TypeError):
            self.L(lcmath.ε)
        with self.assertRaises(TypeError):
            self.L(2j)

    def test_math_cases(self):
        self.assertMathTestCases('math_testcases.txt')

    def test_eps(self):
        self.assertNotEqual(self.ε, 0)
        self.assertNotEqual(self.ε, self.L(0))
        self.assertTrue(self.ε)
        self.assertFalse(self.L(0))
        self.assertFalse(self.L(2)+self.ε-self.L(2)-self.ε)
        self.assertEqual(0, self.L(0))
        self.assertEqual(self.ε+self.ε, 2*self.ε)
        self.assertEqual(self.ε*2, 2*self.ε)
        self.assertEqual(self.L(2)*self.ε, 2*self.ε)
        self.assertEqual(self.ε*self.L(2), 2*self.ε)
        self.assertEqual(self.L(1)+self.ε, self.ε+self.L(1))
        self.assertEqual(self.L(1)+self.ε, self.L(1, 0, ((0, 1), (1, 1))))
        self.assertEqual((self.L(1)+self.ε).st(), 1)
        self.assertEqual(self.ε * (1/self.ε), 1)
        self.assertEqual(self.ε*self.ε, self.ε**2)
        self.assertEqual(self.ε*self.ε, pow(self.ε, 2))
        self.assertEqual(1/(1-self.ε), sum(self.ε**n for n in range(self.ε._TERMS)))
        self.assertEqual(self.exp(self.ε),
                         sum((1/math.factorial(i))*self.ε**i
                             for i in range(self.ε._TERMS)))
        self.assertClose(self.cos(self.sin(self.cos(self.sin(self.ε)))),
                         .666367, rel_tol=1e-6)        
        self.assertEqual((1+self.ε)**self.pi, 1)
        # This is really testing internal implementaton details, so it's
        # probably not a great test...
        self.assertEqual(self.ε**self.pi, self.ε**(2**-48))
        self.assertClose(self.exp(self.log(2+3*self.ε)),
                         self.log(self.exp(2+3*self.ε)))
        with self.assertRaises(OverflowError):
            (1/self.ε).exp()
        with self.assertRaises(OverflowError):
            2**(1/self.ε)

    def test_sign(self):
        self.assertEqual((1-self.ε+self.ε**2).copysign(-1), -1+self.ε-self.ε**2)
        self.assertEqual((-1+self.ε-self.ε**2).copysign(-1), -1+self.ε-self.ε**2)
        self.assertEqual((1-self.ε+self.ε**2).copysign(self.L(-0.0)), -1+self.ε-self.ε**2)
        self.assertEqual(math.copysign(1, self.L(0.0).copysign(-1).front), -1)

    def test_compare(self):
        self.assertLess(self.L(1), self.L(2))
        self.assertLess(self.ε, 1e-300)
        self.assertLess(0, self.ε)
        self.assertLess(-self.ε, self.ε)
        self.assertLess(-1e-300, self.ε)
        self.assertLess(-1e-300, -self.ε)
        self.assertLess(self.ε**2, self.ε)
        self.assertLess(self.ε, self.ε*2)
        self.assertLess(self.ε, self.sqrt(self.ε))
        self.assertLess(sys.float_info.max, 1/self.ε)
        
    def test_st(self):
        self.assertEqual(self.st(self.ε), 0)
        self.assertEqual(self.st(1+2*self.ε), 1)
        self.assertEqual(self.st((1+self.ε)*self.ε), 0)
        
    def test_str(self):
        self.assertEqual(str(self.L(1)+self.ε), "1+ε")
        self.assertEqual(str(1-self.ε), "1-ε")
        self.assertEqual(str(1-self.ε**3), "1-ε**3")
        self.assertEqual(str((1-self.ε)**2), "1-2*ε+ε**2")
        self.assertEqual(str(1 + 1/self.ε), "ε**-1+1")
        self.assertEqual(str(1-self.ε**3-self.ε**-3), "-ε**-3+1-ε**3")
        
    def test_underflow_overflow(self):
        # TODO: handle once we fix underflow
        min = sys.float_info.min
        self.assertNotEqual(min * self.ε, 0)
        #self.assertNotEqual(min * self.ε * min, 0)
        #self.assertNotEqual(min * (min * self.ε), 0)
        #self.assertNotEqual((min + self.ε)**2, 0)

        # TODO: handle once we fix overflow
        max_plus = sys.float_info.max + self.ε
        # currently gives inf+inf*self.ε
        # (inf calculator gives inf + nan/nan*self.ε + nan/nan*self.ε**2)
        #self.assert??((max+self.ε)**2, ???)

class TestLeviCivitaComplex(_TestBaseLeviCivita):
    _MATHMOD = cmath
    _MATHCLASS = complex
    _LCMOD = lcmath
    _LCCLASS = lcmath.LeviCivitaComplex
    
    def assertSmallerMagnitude(self, first, second, msg=None):
        mfirst, msecond = abs(self.L(first)), abs(self.L(second))
        if ((-mfirst.leading, mfirst.front, mfirst.series) <
            (-msecond.leading, msecond.front, msecond.series)):
            return
        sr = unittest.util.safe_repr
        standardMsg = f'|{first}| >= |{second}| (|{sr(mfirst)}| >= |{sr(msecond)}|)'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)
    
    def test_real(self):
        self.assertFalse(self.L(2)+self.L(-2))
        self.assertEqual(self.L(2)+self.L(-3), self.L(-1))
        self.assertEqual(self.L(2)+(-3), self.L(-1))
        self.assertEqual(self.L(2)-3, self.L(-1))
        self.assertEqual(2+self.L(-3), self.L(-1))
        self.assertEqual(2-self.L(3), self.L(-1))
        self.assertEqual(self.L(2)+self.L(-3), -1)
        with self.assertRaises(ZeroDivisionError):
            self.L(2)/self.L(0)

    def test_complex(self):
        self.assertFalse(self.L(2)+self.L(-2))
        self.assertEqual(self.L(2)+self.L(-3), self.L(-1))
        self.assertEqual(self.L(2)+(-3), self.L(-1))
        self.assertEqual(self.L(2)-3, self.L(-1))
        self.assertEqual(2+self.L(-3), self.L(-1))
        self.assertEqual(2-self.L(3), self.L(-1))
        self.assertEqual(self.L(2)+self.L(-3), -1)
        with self.assertRaises(ZeroDivisionError):
            self.L(2)/self.L(0)
        self.assertEqual(self.L(1+2j) * self.L(2+1j), self.L(5j))
        self.assertClose(self.sqrt(self.L(-1)), 1j)
        self.assertClose(self.L(-1)**.5, 1j)
        self.assertClose(self.L(-1)**self.L(.5), 1j)
        with self.assertRaises(ZeroDivisionError):
            (2+3j)/self.L(0j)

    def test_cmath_cases(self):
        self.assertMathTestCases('cmath_testcases.txt')

    def test_eps(self):
        self.assertNotEqual(self.ε, 0)
        self.assertNotEqual(self.ε, self.L(0))
        self.assertTrue(self.ε)
        self.assertFalse(self.L(0))
        self.assertFalse(self.L(2)+self.ε-self.L(2)-self.ε)
        self.assertEqual(0, self.L(0))
        self.assertEqual(self.ε+self.ε, 2*self.ε)
        self.assertEqual(self.ε*2, 2*self.ε)
        self.assertEqual(self.L(2)*self.ε, 2*self.ε)
        self.assertEqual(self.ε*self.L(2), 2*self.ε)
        self.assertEqual(self.L(1)+self.ε, self.ε+self.L(1))
        self.assertEqual(self.L(1)+self.ε, self.L(1, 0, ((0, 1), (1, 1))))
        self.assertEqual((self.L(1)+self.ε).st(), 1)
        self.assertEqual(self.ε * (1/self.ε), 1)
        self.assertClose((1+2j+self.ε) * (2+1j+self.ε), 5j+(3+3j)*self.ε+self.ε**2)
        self.assertEqual(self.ε*self.ε, self.ε**2)
        self.assertEqual(self.ε*self.ε, pow(self.ε, 2))
        self.assertEqual(1/(1-self.ε), sum(self.ε**n for n in range(self.ε._TERMS)))
        self.assertEqual(self.exp(self.ε),
                         sum((1/math.factorial(i))*self.ε**i
                             for i in range(self.ε._TERMS)))
        self.assertClose(self.cos(self.sin(self.cos(self.sin(self.ε)))),
                         .666367, rel_tol=1e-6)        
        self.assertEqual((1+self.ε)**self.pi, 1)
        # This is really testing internal implementaton details, so it's
        # probably not a great test...
        self.assertEqual(self.ε**self.pi, self.ε**(2**-48))
        self.assertClose(self.exp(self.log(2+3*self.ε)),
                         self.log(self.exp(2+3*self.ε)))
        with self.assertRaises(OverflowError):
            (1/self.ε).exp()
        with self.assertRaises(OverflowError):
            2**(1/self.ε)

    def test_compare(self):
        self.assertSmallerMagnitude(self.L(1), self.L(2))
        self.assertSmallerMagnitude(self.ε, 1e-300)
        self.assertSmallerMagnitude(self.ε, -1e-300)
        self.assertSmallerMagnitude(-self.ε, 1e-300)
        self.assertSmallerMagnitude(self.ε**2, self.ε)
        self.assertSmallerMagnitude(self.ε, self.ε*2)
        self.assertSmallerMagnitude(self.ε, self.sqrt(self.ε))
        
    def test_st(self):
        self.assertEqual(self.st(self.ε), 0)
        self.assertEqual(self.st(1+2j), 1+2j)
        self.assertEqual(self.st(self.L(1+2j)), 1+2j)
        self.assertEqual(self.st(1+2j+self.ε), 1+2j)
        self.assertEqual(self.st(1+2j+self.ε*2j), 1+2j)
        self.assertEqual(self.st((1+self.ε)*self.ε), 0)
        self.assertEqual(self.st((1+2j+self.ε) * (2+1j+self.ε)), 5j)
        
    def test_str(self):
        self.assertEqual(str(self.L(1)+self.ε), "(1+0j)+ε")
        self.assertEqual(str(1-self.ε), "(1+0j)-ε")
        self.assertEqual(str(1-self.ε**3), "(1+0j)-ε**3")
        self.assertEqual(str((1-self.ε)**2), "(1+0j)-(2+0j)*ε+ε**2")
        self.assertEqual(str(1 + 1/self.ε), "ε**-1+(1+0j)")
        self.assertEqual(str(1-self.ε**3-self.ε**-3), "-ε**-3+(1+0j)-ε**3")
        
    # TODO add order comparisons once reals are added
    def test_underflow_overflow(self):
        # TODO: handle once we fix underflow
        min = sys.float_info.min
        self.assertNotEqual(min * self.ε, 0)
        #self.assertNotEqual(min * self.ε * min, 0)
        #self.assertNotEqual(min * (min * self.ε), 0)
        #self.assertNotEqual((min + self.ε)**2, 0)

        # TODO: handle once we fix overflow
        max_plus = sys.float_info.max + self.ε
        # currently gives inf+inf*self.ε
        # (inf calculator gives inf + nan/nan*self.ε + nan/nan*self.ε**2)
        #self.assert??((max+self.ε)**2, ???)
        
if __name__ == '__main__':
    unittest.main()
