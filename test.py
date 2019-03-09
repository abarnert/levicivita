import cmath
import math
import sys
import unittest

from levicivita import cpx as c
from levicivita import real as r

C = c.LeviCivitaComplex
ceps = c.eps
F = r.LeviCivitaFloat
feps = r.eps

class TestLeviCivitaFloat(unittest.TestCase):
    def assertClose(self, first, second, msg=None, *,
                    rel_tol=1e-09, abs_tol=0.0):
        if r.isclose(first, second, rel_tol=rel_tol, abs_tol=abs_tol):
            return
        diff = abs(first - second)
        sr = unittest.util.safe_repr
        stdiff = diff if r.isinf(diff) else r.st(diff)
        standardMsg = f'{first} != {second} ({sr(first)} != {sr(second)}) within rel_tol {rel_tol}, abs_tol {abs_tol} ({stdiff} difference)'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)

    def assertLess(self, first, second, msg=None):
        if first < second:
            return
        sr = unittest.util.safe_repr
        standardMsg = f'{first} >= {second} ({sr(first)} >= {sr(second)})'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)
    
    def test_real(self):
        self.assertFalse(F(2)+F(-2))
        self.assertEqual(F(2)+F(-3), F(-1))
        self.assertEqual(F(2)+(-3), F(-1))
        self.assertEqual(F(2)-3, F(-1))
        self.assertEqual(2+F(-3), F(-1))
        self.assertEqual(2-F(3), F(-1))
        self.assertEqual(F(2)+F(-3), -1)
        with self.assertRaises(ZeroDivisionError):
            F(2)/F(0)

    def test_complex(self):
        with self.assertRaises(TypeError):
            F(C(2))
        with self.assertRaises(TypeError):
            F(ceps)
        with self.assertRaises(TypeError):
            F(2j)

    def test_eps(self):
        self.assertNotEqual(feps, 0)
        self.assertNotEqual(feps, F(0))
        self.assertTrue(feps)
        self.assertFalse(F(0))
        self.assertFalse(F(2)+feps-F(2)-feps)
        self.assertEqual(0, F(0))
        self.assertEqual(feps+feps, 2*feps)
        self.assertEqual(feps*2, 2*feps)
        self.assertEqual(F(2)*feps, 2*feps)
        self.assertEqual(feps*F(2), 2*feps)
        self.assertEqual(F(1)+feps, feps+F(1))
        self.assertEqual(F(1)+feps, F(1, 0, ((0, 1), (1, 1))))
        self.assertEqual((F(1)+feps).st(), 1)
        self.assertEqual(feps * (1/feps), 1)
        self.assertEqual(feps*feps, feps**2)
        self.assertEqual(feps*feps, pow(feps, 2))
        self.assertEqual(1/(1-feps), sum(feps**n for n in range(feps._TERMS)))
        self.assertEqual(r.exp(feps),
                         sum((1/math.factorial(i))*feps**i
                             for i in range(feps._TERMS)))
        self.assertClose(r.cos(r.sin(r.cos(r.sin(feps)))), .666367,
                         rel_tol=1e-6)
        self.assertEqual((1+feps)**r.pi, 1)
        # This is really testing internal implementaton details, so it's
        # probably not a great test...
        self.assertEqual(feps**r.pi, feps**(2**-48))
        self.assertClose(r.exp(r.log(2+3*feps)), r.log(r.exp(2+3*feps)),
                         rel_tol=1e-6)
        with self.assertRaises(OverflowError):
            2**(1/feps)

    def test_sign(self):
        self.assertEqual((1-feps+feps**2).copysign(-1), -1+feps-feps**2)
        self.assertEqual((-1+feps-feps**2).copysign(-1), -1+feps-feps**2)
        self.assertEqual((1-feps+feps**2).copysign(F(-0.0)), -1+feps-feps**2)
        self.assertEqual(math.copysign(1, F(0.0).copysign(-1).front), -1)

    def test_compare(self):
        self.assertLess(F(1), F(2))
        self.assertLess(feps, 1e-300)
        self.assertLess(0, feps)
        self.assertLess(-feps, feps)
        self.assertLess(-1e-300, feps)
        self.assertLess(-1e-300, -feps)
        self.assertLess(feps**2, feps)
        self.assertLess(feps, feps*2)
        self.assertLess(feps, r.sqrt(feps))
        self.assertLess(sys.float_info.max, 1/feps)
        
    def test_st(self):
        self.assertEqual(r.st(feps), 0)
        self.assertEqual(r.st(1+2*feps), 1)
        self.assertEqual(r.st((1+feps)*feps), 0)
        
    def test_str(self):
        self.assertEqual(str(F(1)+feps), "1+ε")
        self.assertEqual(str(1-feps), "1-ε")
        self.assertEqual(str(1-feps**3), "1-ε**3")
        self.assertEqual(str((1-feps)**2), "1-2*ε+ε**2")
        self.assertEqual(str(1 + 1/feps), "ε**-1+1")
        self.assertEqual(str(1-feps**3-feps**-3), "-ε**-3+1-ε**3")
        
    def test_underflow_overflow(self):
        # TODO: handle once we fix underflow
        min = sys.float_info.min
        self.assertNotEqual(min * feps, 0)
        #self.assertNotEqual(min * feps * min, 0)
        #self.assertNotEqual(min * (min * feps), 0)
        #self.assertNotEqual((min + feps)**2, 0)

        # TODO: handle once we fix overflow
        max_feps = sys.float_info.max + feps
        # currently gives inf+inf*feps
        # (inf calculator gives inf + nan/nan*feps + nan/nan*feps**2)
        #self.assert??((max+feps)**2, ???)

class TestLeviCivitaComplex(unittest.TestCase):
    def assertClose(self, first, second, msg=None, *,
                    rel_tol=1e-09, abs_tol=0.0):
        if c.isclose(first, second, rel_tol=rel_tol, abs_tol=abs_tol):
            return
        diff = abs(first - second)
        sr = unittest.util.safe_repr
        stdiff = diff if c.isinf(diff) else st(diff)
        standardMsg = f'{first} != {second} ({sr(first)} != {sr(second)}) within rel_tol {rel_tol}, abs_tol {abs_tol} ({stdiff} difference)'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)

    def assertSmallerMagnitude(self, first, second, msg=None):
        mfirst, msecond = abs(C(first)), abs(C(second))
        if ((-mfirst.leading, mfirst.front, mfirst.series) <
            (-msecond.leading, msecond.front, msecond.series)):
            return
        sr = unittest.util.safe_repr
        standardMsg = f'|{first}| >= |{second}| (|{sr(mfirst)}| >= |{sr(msecond)}|)'
        msg = self._formatMessage(msg, standardMsg)
        raise self.failureException(msg)
    
    def test_real(self):
        self.assertFalse(C(2)+C(-2))
        self.assertEqual(C(2)+C(-3), C(-1))
        self.assertEqual(C(2)+(-3), C(-1))
        self.assertEqual(C(2)-3, C(-1))
        self.assertEqual(2+C(-3), C(-1))
        self.assertEqual(2-C(3), C(-1))
        self.assertEqual(C(2)+C(-3), -1)
        with self.assertRaises(ZeroDivisionError):
            C(2)/C(0)

    def test_complex(self):
        self.assertFalse(C(2)+C(-2))
        self.assertEqual(C(2)+C(-3), C(-1))
        self.assertEqual(C(2)+(-3), C(-1))
        self.assertEqual(C(2)-3, C(-1))
        self.assertEqual(2+C(-3), C(-1))
        self.assertEqual(2-C(3), C(-1))
        self.assertEqual(C(2)+C(-3), -1)
        with self.assertRaises(ZeroDivisionError):
            C(2)/C(0)
        self.assertEqual(C(1+2j) * C(2+1j), C(5j))
        self.assertClose(c.sqrt(C(-1)), 1j)
        self.assertClose(C(-1)**.5, 1j, rel_tol=1e-06)
        self.assertClose(C(-1)**C(.5), 1j, rel_tol=1e-06)
        with self.assertRaises(ZeroDivisionError):
            (2+3j)/C(0j)

    def test_eps(self):
        self.assertNotEqual(ceps, 0)
        self.assertNotEqual(ceps, C(0))
        self.assertTrue(ceps)
        self.assertFalse(C(0))
        self.assertFalse(C(2)+ceps-C(2)-ceps)
        self.assertEqual(0, C(0))
        self.assertEqual(ceps+ceps, 2*ceps)
        self.assertEqual(ceps*2, 2*ceps)
        self.assertEqual(C(2)*ceps, 2*ceps)
        self.assertEqual(ceps*C(2), 2*ceps)
        self.assertEqual(C(1)+ceps, ceps+C(1))
        self.assertEqual(C(1)+ceps, C(1, 0, ((0, 1), (1, 1))))
        self.assertEqual((C(1)+ceps).st(), 1)
        self.assertEqual(ceps * (1/ceps), 1)
        self.assertClose((1+2j+ceps) * (2+1j+ceps), 5j+(3+3j)*ceps+ceps**2)
        self.assertEqual(ceps*ceps, ceps**2)
        self.assertEqual(ceps*ceps, pow(ceps, 2))
        self.assertEqual(1/(1-ceps), sum(ceps**n for n in range(ceps._TERMS)))
        self.assertEqual(c.exp(ceps),
                         sum((1/math.factorial(i))*ceps**i
                             for i in range(ceps._TERMS)))
        self.assertClose(c.cos(c.sin(c.cos(c.sin(ceps)))), .666367,
                         rel_tol=1e-6)        
        self.assertEqual((1+ceps)**c.pi, 1)
        # This is really testing internal implementaton details, so it's
        # probably not a great test...
        self.assertEqual(ceps**c.pi, ceps**(2**-48))
        self.assertClose(c.exp(c.log(2+3*ceps)), c.log(c.exp(2+3*ceps)),
                         rel_tol=1e-6)
        with self.assertRaises(OverflowError):
            2**(1/ceps)

    def test_compare(self):
        self.assertSmallerMagnitude(C(1), C(2))
        self.assertSmallerMagnitude(ceps, 1e-300)
        self.assertSmallerMagnitude(ceps, -1e-300)
        self.assertSmallerMagnitude(-ceps, 1e-300)
        self.assertSmallerMagnitude(ceps**2, ceps)
        self.assertSmallerMagnitude(ceps, ceps*2)
        self.assertSmallerMagnitude(ceps, c.sqrt(ceps))
        
    def test_st(self):
        self.assertEqual(c.st(ceps), 0)
        self.assertEqual(c.st(1+2j), 1+2j)
        self.assertEqual(c.st(C(1+2j)), 1+2j)
        self.assertEqual(c.st(1+2j+ceps), 1+2j)
        self.assertEqual(c.st(1+2j+ceps*2j), 1+2j)
        self.assertEqual(c.st((1+ceps)*ceps), 0)
        self.assertEqual(c.st((1+2j+ceps) * (2+1j+ceps)), 5j)
        
    def test_str(self):
        self.assertEqual(str(C(1)+ceps), "1+ε")
        self.assertEqual(str(1-ceps), "1-ε")
        self.assertEqual(str(1-ceps**3), "1-ε**3")
        self.assertEqual(str((1-ceps)**2), "1-2*ε+ε**2")
        self.assertEqual(str(1 + 1/ceps), "ε**-1+1")
        self.assertEqual(str(1-ceps**3-ceps**-3), "-ε**-3+1-ε**3")
        
    # TODO add order comparisons once reals are added
    def test_underflow_overflow(self):
        # TODO: handle once we fix underflow
        min = sys.float_info.min
        self.assertNotEqual(min * ceps, 0)
        #self.assertNotEqual(min * ceps * min, 0)
        #self.assertNotEqual(min * (min * ceps), 0)
        #self.assertNotEqual((min + ceps)**2, 0)

        # TODO: handle once we fix overflow
        max_ceps = sys.float_info.max + ceps
        # currently gives inf+inf*ceps
        # (inf calculator gives inf + nan/nan*ceps + nan/nan*ceps**2)
        #self.assert??((max+ceps)**2, ???)
        
if __name__ == '__main__':
    unittest.main()
