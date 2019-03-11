#!/usr/bin/env python3

import code
from textwrap import dedent

console = code.InteractiveConsole()
for line in dedent('''\
    import readline
    import levicivita
    from levicivita.real import *
    from levicivita import cpx
    from levicivita import real
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
    '''))
