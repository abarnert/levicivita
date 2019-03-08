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

The main user interface is:

 * constant `ε`
  * also named named `eps` and `epsilon` (and `d`, for compatibility
    with the `inf` calculator)
  * `1/ε`, `2*ε`, `ε**2`, `1+ε`, `1+2j+(1-2j)*ε`, etc. give you any
    other L-C numbers you need
 * all `cmath`-style functions, extended to deal with L-C numbers 
  * `phase(x)`, `polar(x)`, `rect(r, phi)`
  * `exp(x)`, `log(x[, base])`, `log10(x)`, `sqrt(x)`
  * `acos(x)`, `asin(x)`, `atan(x)`, `cos(x)`, `sin(x)`, `tan(x)`
  * `acosh(x)`, `asinh(x)`, `atanh(x)`, `cosh(x)`, `sinh(x)`, `tanh(x)`
  * `isfinite(x)`, `isinf(x)`, `isnan(x)`
  * `isclose(a, b, *, rel_tol=1-09, abs_tol=0.0)`
 * `st(x)` returns the standard (non-infinitestimal) part of `x`,
   raising an exception on infinite `x`
 * `change_terms(n)` to change the approximation depth

Most of these should be obvious, except the last two.

Two infinite numbers are close if they have the same leading exponent
with coefficients that are close (ignoring `abs_tol`). Two finite
numbers are close if their standard parts are close. (Note that this
means any two infinitesimals are always close.)

L-C numbers are approximated to the first N terms (and of course to
IEEE double values for coefficients), by default 5. This can be
changed globally, by calling `change_terms` with a different `N`.

You almost always want to print out the `str` rather than `repr` of
these numbers.

You should rarely need to construct `LeviCivitaComplex` numbers
directly, but if you need to:

 * `LeviCivitaComplex()` == `0`
 * `LeviCivitaComplex(c)` == `c` for any `complex` or compatible
   number, or any existing `LeviCivitaComplex` number.
 * `LeviCivitaComplex(front, leading)` == `front * ε**leading`.
 * `LeviCivitaComplex(front, leading, series)`: see **Implementation**.

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

Currently, the package only handles complex L-C numbers, not
real. Which means that it doesn't handle ordering at all (even though
the simple ordering properties are one of the reasons L-C is so nice).

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
do that.)

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

# Implementation

The implementation was borrowed from the JavaScript browser calculator
[`inf`](http://www.lightandmatter.com/calc/inf/) (available
[on GitHub](https://github.com/bcrowell/inf)). However, it's a new
implementation of the same structures--on the plus side, that means
it's a lot more Pythonic than a straight port would have been; on the
minus side, that means it doesn't handle all of the features that a
port would have.

Internally, L-C numbers are implemented by the class
`LeviCivitaComplex`, which normalizes numbers as:

 * `front`: standard part (`complex`)
 * `leading`: leading epsilon exponent (`int` or `Fraction`)
 * `series`: sequence of up to `TERMS` pairs of `(a, q)` values
   (`real` or `complex`, `int` or `Fraction`) for additional epsilon
   terms

This idea was borrowed straight from `inf`. It makes a lot of things
(like comparing standard parts) very easy, at the cost of making it
more complicated to construct numbers. Numbers are always normalized
on construction. (I believe that's not true of `inf`; you can
construct a number any way you want, and calling a `tidy` method is
optional.) 

As with other numeric types, `LeviCivitaComplex` values are immutable
(which also isn't true for `inf`), and for two equal numbers, neither
identity nor nonidentity is guaranteed. (In other words, `__new__` may
return an existing number if it's convenient.)

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

The class fully implements `numbers.Complex`, with all of the
operators doing exactly what you'd expect (including the reverse
forms). It plays nicely with any type in the numeric tower--including
the builtin types and `fractions.Fraction`, but not `decimal.Decimal`.

It also has methods for all of the module functions in `cmath`; the
module functions call those methods when given `LeviCivitaComplex`
objects and defer to `cmath` otherwise.

`1/x`, `exp(x)`, and `log(x)` are implemented through Taylor series,
which are statically built (to `TERMS` terms) when you call
`change_terms`. `pow(x, y)` is implemented through a Taylor series
built on the fly (although it is optimized for integral `y`). All of
the other `cmath` functions are built in terms of these. For example,
`acos(x)` is `-j * log(x + sqrt(1-x**2)`. Note that these
implementations are not always as precise (not to mention efficient)
as the ones in `cmath` and `math`.

# History

 * 0.0.1 2019-03-08: initial implementation

# TODO

 * Handle `float` as well as `complex` (including adding comparisons).
 * `conjugate`, `real`, `imag`, `phase`, `rect`, `polar` are missing.
 * Extend `isclose` with `eps_tol` or something to handle terms up to
   `eps**n`? (It would certainly be useful for unit testing, if
   nothing else. For example, `log(exp(2+3*ε)) == exp(log(2+3*ε))` is
   probably false; `isclose` only tests that the `2` is correct; we
   really want to test that the `2` is close to `2`, the `3` is close
   to `3`, and the next few terms are all close to `0`.
 * Are `exp(1/ε)`, `log(1/ε)`, `log(ε)`, etc. really out of range? If
   not, implement them.
 * Should `isinf` and friends treat infinite numbers like `inf`? And
   should there even be an `inf = cmath.inf` instead of, say, `1/ε`?
 * Replace `change_terms` with something more like `decimal` contexts.
 * Maybe `__float__` and `__int__` should raise custom `TypeError`
   like `complex` does?
 * Separate out unit tests.
 * Add more tests and more error handling.
 * Can I look deeper into implementation of `inf` despite GPL? There
   are a lot of places where I guessed at the implementation from the
   function names, and may have gotten things wrong...
 * Look at construction from `Fraction`. (Disallow? Convert to `float`?)
 * Consider working with `Decimal`.
 * Construct from string.
 * Look over handling of `inf`, `nan`, etc., especially
   overflow/underflow.
 * Improve normalization: Don't sort `series` multiple times (maybe
   don't sort at all; try collecting via `dict`/`Counter`); don't
   convert back to `tuple` multiple times; etc.
 * Consider implementing `pow` in terms of `exp` and `log`.
 * Handle cases where `pow` (and maybe `exp`, etc.?) gives complex for
   complex-typed but real-valued arguments, ideally as well as `cmath`
   does.
 * Compare errors for different implementations of transcendental
   functions (and maybe even use explicit Taylor series for each instead
   of using `exp`).
