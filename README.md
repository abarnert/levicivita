# levicivita
Levi-Civita infinite and infinitesimal values

# Introduction

A Levi-Civita field extends the real or complex numbers with infinite
and infinitesimal numbers.

They're typically used for automatic differentiation, rather than
directly calculating with them. But it can be instructive to play with
infinitesimal and infinite numbers, hence this package.

Any Levi-Civita number is a sum of `aε^q` values, where each a is real
or complex, and each q is rational. Basic arithmetic works in the
obvious way, and ordering is just lexicographical (dictionary)
ordering (if you keep the terms sorted in ascending `q`).

 * A real number `r` is just `r*(ε**0)`
 * A complex number `x+yj` is just `(x+yj)*(ε**0)`
 * An infinitesimal looks like `ε` or `7*ε**(3/2)` (or `0`, of course)
 * An infinite number looks like `1/ε` or `ε**-2 + 2*ε**-1 + 1`
 * A number differing infinitesimally from a real looks like `7+ε`, or
   maybe `7 + 49*ε + 14*(ε**2)`

Infinitesimals built from `ε` and infinite numbers built from `1/ε`
work exactly the way you'd expect:

 * `0 < ε < MIN_FLOAT` (smaller than any positive real number)
 * `1-ε < 1 < 1+ε`
 * `ε/2 < ε < 7*ε`
 * `ε**2 < ε < sqrt(ε)`
 * `MAX_FLOAT < 1/ε < INF` (larger than any real number, but not
   "absolute infinity")
 * `sqrt(1/ε) < 1/ε < (1/ε)**2`

For more detail, see **Math** below.

# Usage

Interactively, either:

 * `from levicivita.real import *` or
 * `from levicivita.cpx import *`

For non-interactive use, you probably don't want `import *` of course.

The main user interface is:

 * constant `ε`
  * also named named `eps` and `epsilon` (and `d`, for compatibility
    with the `inf` calculator)
  * `1/ε`, `2*ε`, `ε**2`, `1+ε`, `1+2j+(1-2j)*ε`, etc. give you any
    other L-C numbers you need
 * `math`/`cmath`-style functions, extended to deal with L-C numbers 
  * `phase(x)`, `polar(x)`, `rect(r, phi)` (complex only)
  * `ceil(x)`, `floor(x)`, `trunc(x)`, `round(x, digits)` (real only)
  * `exp(x)`, `log(x[, base])`, `log10(x)`, `sqrt(x)`
  * `acos(x)`, `asin(x)`, `atan(x)`, `cos(x)`, `sin(x)`, `tan(x)`
  * `acosh(x)`, `asinh(x)`, `atanh(x)`, `cosh(x)`, `sinh(x)`, `tanh(x)`
  * `isfinite(x)`, `isinf(x)`, `isnan(x)`
  * `isclose(a, b, *, rel_tol=1-09, abs_tol=0.0)`
 * `st(x)` returns the standard (non-infinitestimal) part of `x`,
   raising an exception on infinite `x`
 * `change_terms(n)` to change the approximation depth

Most of these should be obvious, except `isclose` and `change_terms`.

Two infinite numbers are close if they have the same leading exponent
with coefficients that are close (ignoring `abs_tol`). Two finite
numbers are close if their standard parts are close (by both `rel_tol`
and `abs_tol`). Note that this means that two infinitesimals are
always close.

L-C numbers are approximated to the first N terms (and of course to
IEEE double values for coefficients), by default 5 (although twice as
many are used for intermediate calculations). This can be changed
globally, by calling `change_terms` with a different `N`.

You almost always want to print out the `str` rather than `repr` of
these numbers.

You should rarely need to construct Levi-Civita numbers directly, but
if you need to:

 * `LeviCivitaComplex()` == `0`
 * `LeviCivitaComplex(c)` == `c` for any `complex` or compatible
   number, or any existing `LeviCivitaComplex` number.
 * `LeviCivitaComplex(front, leading)` == `front * ε**leading`.
 * `LeviCivitaComplex(front, leading, series)`: see **Implementation**.

... and likewise for `LeviCivitaFloat` (except of course that the
`front` and the `series` coefficients must be real).

# Automatic differentiation

Given any Python function, you can calculate its derivative at any
real or complex value. After all, the derivative is defined as
`(f(x+ε) - f(x)) / ε`, and we can directly approximate that
calculation. (In fact, this is the main use for Levi-Civita
infinitesimals.)

To play with this:

    import diff
	def func(x):
	    return 2*x**3 + 3*x**2 + 4*x + 5
	print(diff.derivative(func, 3))
	print(diff.derivative(func, 3j))
	
This should print out something very close to `76` and `-50+18j`, the
same values you'd get by manually implementing the derivative function
as `6*x**2 + 6*x + 4`.

However, there is a big caveat. The function above works because the
only thing it does with `x` is apply operators to it. If you tried to,
say, call `math.gamma`, or any other function (that isn't recursively
ultimately defined in terms of operators), you'd just get a
`TypeError`. Sometimes you can remove the function call (e.g., replace
`math.exp(x)` with `math.e ** x`, or `cmath.sqrt(x)` with `x**.5`). If
the function you want to call is one of the ones provided by the
`levicivita.real` package, you can use that in place of `math` (and
likewise for `cmath`); you can even monkeypatch some existing code's
globals to do that if you really want to. But if it's calling some
other function (or some C extension like `numpy`), there's really no
way around that.

Of course if you really want to do automatic differentiation, you
could use a C implementation that hooks `libmath`. I believe
[`COSY`](https://www.bt.pa.msu.edu/index_cosy.htm) offers that
functionality, for example.

But this should be enough to play with the idea.

# Math

Levi-Civita infinitesimals are not as powerful as hyperreals, but
they're a lot more useful for computation.

They can't extend arbitrary real functions, but they can extend
functions with Taylor series expansions, which is often good enough.

The advantage of L-C is that you can approximate numbers and functions
very nicely by keeping just the first few terms, with approximate
coefficients (say, IEEE doubles). Mathematically, it's guaranteed that
the support (the set of indices where `a` is nonzero) is a left-finite
set. It's _not_ guaranteed that the coefficients grow fast enough that
you can ignore all but the first few, or that you're not going to do
something that cancels out the first 350 so everything depends on the
351st, and you can create pathological cases that break the
approximation. But in practice, you shouldn't need to worry. In
general, you can just use them with the same caveats you worry about
for `float`.

For example, approximating the `exp` function gives useful values with
only the first few terms and only IEEE coefficients:

 * `e^ε = 1 + ε + ε^2/2 + ... + ε^n/n!`

The 10th term and beyond have IEEE rounding errors, but unless you do
something that cancels out the first 9 terms, that won't even show up,
and unless you do something that cancels out the first 9 terms and
needs 17 digits of precision, it won't be a problem.

See [Wikipedia](https://en.wikipedia.org/wiki/Levi-Civita_field) and
["Analysis on the Levi-Civita Field: A Brief Overview"](http://www2.physics.umanitoba.ca/u/khodr/Publications/RS-Overview-offprints.pdf)
by Shamseddine and Berz for more information.

# Design Limitations

The implementation of transcendental methods on L-C float is based
heavily on complex. For `exp`/`log`/`pow` we do the complex math,
check that the imaginary part is small enough to be a plausible error
(`1e-15` in the `front`), and throw it away, which is pretty
hacky. For `sin` and friends, we don't even do the check (because
in practice, `sin(cos(sin(ε)))` would already have much more than
`1e-15` error), which is even hackier. There is probably a better
solution.

Using complex coefficients of each term, instead of having separate
real and imaginary series, means a few things that IEEE complex
numbers handle nicely (like numbers with finite real but infinite
imaginary parts) aren't as nice. I think the benefits outweigh the
costs, however.

It needs better support for dealing with IEEE infinity and NaN in
general, and for overflow and underflow in particular.

The internal structure (which is based on the documentation of the
`inf` calculator; see below) may not be ideal. For example, do we
really need unbounded exact fractions for the exponents? Probably not,
and it probably slows things down. (I suspect `inf` doesn't actually
do that, but I could be wrong.)

You almost always want to print out the `str` rather than the `repr`,
which is annoying for interactive use. Maybe there should be a way to
switch this without monkeypatching, maybe even always on by default in
interactive mode?

The number of series terms is a global (class) constant, which can
only be set by calling `change_terms` on the class. This should
probably be more like a `decimal` module context.

The class doesn't interoperate very nicely with `decimal`. In part,
this is because it's not clear how it _should_. Would we want to keep
all of the coefficients as decimals with the specified precision? What
about the epsilon exponents?

The design of `isclose` tries maybe too hard to match `cmath`. As
written, it means any two infinitesimals are always close to each
other. Maybe `rel_tol` should always be applied to `front`, instead of
only when nonpositive. And, even beyond that, it might be useful (at
least for testing this module!) to have a way to specify closeness in
powers (`eps_tol`?) of `ε`, so, e.g., `1+ε` and `1+7*ε` would fail if
you asked for 2 or more powers.

# Implementation

The implementation was borrowed from the JavaScript browser calculator
[`inf`](http://www.lightandmatter.com/calc/inf/) (available
[on GitHub](https://github.com/bcrowell/inf)). However, it's a new
implementation of the same structures--on the plus side, that means
it's a lot more Pythonic than a straight port would have been; on the
minus side, that means it doesn't handle all of the features that a
port would have.

Internally, L-C numbers are implemented by the classes
`LeviCivitaComplex` and `LeviCivitaFloat`, which normalize numbers as:

 * `front`: standard part (`complex` or `float`, as appropriate)
 * `leading`: leading epsilon exponent (can be `int` or `Fraction`)
 * `series`: sequence of up to `TERMS` pairs of `(a, q)` values
   (`real` or `complex`, `int` or `Fraction`) for additional epsilon
   terms

This idea was borrowed straight from `inf`. It makes a lot of things
(like comparing standard parts) very easy, at the cost of making it
more complicated to construct numbers. Numbers are always normalized
on construction. (I believe that's not true of `inf`; you can
construct a number any way you want, and calling a `tidy` method is
optional.)

As with other numeric types, the values are immutable (which also
isn't true for `inf`), and for two equal numbers, neither identity nor
nonidentity is guaranteed. (In other words, `__new__` may return an
existing number if it's convenient.)

Class attributes `TERMS` and `DISPLAY_TERMS` control how many terms
are preserved at normalization and in `str`, respectively. The
`change_terms(n)` function sets `DISPLAY_TERMS` to `n` and `TERMS` to
`2*n`, and then re-generates the static taylor expansions (see below).

The `repr` will always be a constructor call, which should be
guaranteed to `eval` to the exact same number.

The `str` will always be an expression in coefficients and powers of
`ε`, up to `DISPLAY_TERMS` terms. It will sometimes be `eval`-able to
the same number, but it may be only approximate (e.g., if there are
too many terms, or if there are IEEE rounding issues--and notice that
normalization can cause rounding). It's also less efficient to build
them this way.

The classes fully implement `numbers.Complex` and `numbers.Real`, with
all of the operators doing exactly what you'd expect (including the
reverse forms). They play nicely with any type in the numeric
tower--including the builtin types and `fractions.Fraction`, but not
`decimal.Decimal`.

There are also methods for all of the module functions in `cmath` or
most of the functions in `math`, respectively. The module functions
call those methods when given L-C objects, and defer to `cmath`/`math`
otherwise. (In particular, floats support everything in `math` that's
also in `cmath`, plus `ceil`/`floor`/`round`/`trunc`, but nothing
else, like `degrees` or `erfc`.)

`1/x`, `exp(x)`, and `log(x)` are implemented through Taylor series,
which are statically built (to `TERMS` terms) when you call
`change_terms`. `pow(x, y)` is implemented through a Taylor series
built on the fly (although it is optimized for integral `y`). All of
the other `cmath` functions are built in terms of these. For example,
`acos(x)` is `-j * log(x + sqrt(1-x**2)`. Note that these
implementations are not always as precise (not to mention efficient)
as the ones in `cmath`.

For floats, all of these methods are implemented in terms of the
complex implementation, which is even less precise and efficient
compared to `math`.

# History

 * 0.0.3 2019-03-08: add diff
 * 0.0.2 2019-03-08: add floats, pull tests out
 * 0.0.1 2019-03-08: initial implementation

# BUGS

 * There's clearly something wrong with the transcendental functions.
   For example, `sin(pi/8) - sin(pi/8-d)` should be infinitestimal,
   but it has an error (in the standard part) of `7.8e-8`, which is
   way too big to just discard. This breaks automatic differentiation,
   as it makes the derivative infinite instead of `sqrt(2)`. The good
   news is that the `inf` calculator gives almost the exact same
   results, which implies this may be good enough for playing around
   (as long as you don't try anything like automatic differentiation
   of the transcendental functions). The bad news is that the `inf`
   calculator gives almost the exact same result, which implies that
   there's no better implementation to be found there....

# TODO

 * `conjugate`, `real`, `imag`, `phase`, `rect`, `polar` are missing.
 * `__floordiv__` and `__mod__` are missing.
 * Consider adding other functions from `math`.
 * Is there a better way to implement real transcendental functions?
 * If not, is there at least a better way to refactor things?
 * Consider changing `isclose`: maybe `rel_tol` should apply to
   `front` no matter what, and maybe we want an `eps_tol` or something
   to specify how many powers of epsilon we want to require to be close.
 * Are `exp(1/ε)`, `log(1/ε)`, `log(ε)`, etc. really out of range? If
   not, implement them.
 * Should `isinf` and friends treat infinite numbers like `inf`? And
   should there even be an `inf = cmath.inf` instead of, say, `1/ε`?
 * Replace `change_terms` with something more like `decimal` contexts.
 * Maybe `__float__` and `__int__` should raise custom `TypeError`
   like `complex` does?
 * Add more tests and more error handling.
 * Can I look deeper into implementation of `inf` despite GPL? There
   are a lot of places where I guessed at the implementation from the
   function names, and may have gotten things wrong...
 * Look at construction from `Fraction`. (Disallow? Convert to `float`?)
 * Consider working with `Decimal`.
 * Construct from string?
 * Look over handling of `inf`, `nan`, etc., especially
   overflow/underflow.
 * Improve normalization: Don't sort `series` multiple times (maybe
   don't sort at all; try collecting via `dict`/`Counter`); don't
   convert back to `tuple` multiple times; etc.
 * Handle cases where `pow` (and maybe `exp`, etc.?) gives complex for
   complex-typed but real-valued arguments, ideally as well as `cmath`
   does. Would this help the `_realify` for float (if there's no
   better answer)?
 * Compare errors for different implementations of transcendental
   functions (and maybe even use explicit Taylor series for each instead
   of using `exp`).
