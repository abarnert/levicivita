#!/usr/bin/env python3

import code
import sys
from textwrap import dedent

# TODO: Do we need argparse here?
imported_module = 'lmath'
if len(sys.argv) > 1 and sys.argv[1] == '-c':
    imported_module = 'lcmath'

console = code.InteractiveConsole()
for line in dedent(f'''\
    import readline
    import levicivita
    from levicivita import *
    from levicivita.{imported_module} import *
    from levicivita import lmath
    from levicivita import lcmath
    from levicivita.diff import derivative
    levicivita.LeviCivitaBase._repr = levicivita.LeviCivitaBase.__repr__
    levicivita.LeviCivitaBase.__repr__ = levicivita.LeviCivitaBase.__str__
    ''').splitlines():
    console.push(line)
console.interact(banner=dedent('''\
    This calculator is a full Python interpreter, but with Levi-Civita
    numbers for infinite and infinitesimal values.
    
    The basic infiniesimal is named `ε` or `eps`. Other infinitesimals
    can be created as, e.g., `7*ε` or `ε**2`, and infinite values as
    `1/ε` or `23/ε**2` or `ε**-3`.

    All of the functions from `cmath`, and most of the functions from
    `math`, are available, and handle these numbers appropriately (by
    approximating them to up to 5 powers of ε). Constants from `math`
    like `pi` are also available.

    For example:

        >>> sin(pi/2+ε)
        1+(2.47373e-05)*ε-0.499922*ε**2+0.000149087*ε**3+0.0418552*ε**4
        >>> sin(pi/2+(pi/2)*1j+eps)
        (2.50907+0.000118719j)+(3.01328e-05-2.30049j)*ε-(1.25325+0.00122417j)*ε**2+(0.00244634+0.383416j)*ε**3+(0.105998-0.0016626j)*ε**4
    '''))

