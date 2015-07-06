#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# pyOrbfit extension module
_pyOrbfit = Extension("_pyOrbfit",
                   ["pyOrbfit.i", "fit_radec.c", "orbfit1.c", "nrutil.c", "ephem_earth.c", "aeiderivs.c", "gasdev.c", "abg_to_xyz.c", "gaussj.c",
                    "orbfit2.c", "mrqmin_orbit.c", "abg_to_aei.c", "ludcmp.c", "dms.c", "covsrt.c", "ran1.c",  "lubksb.c", "transforms.c", "mrqcof_orbit.c"],
                   include_dirs = [numpy_include],
                   )


# NumyTypemapTests setup
setup(  name        = "pyOrbfit module",
        description = "Python wrapper for Gary Bernstein's orbfit code",
        author      = "David Gerdes",
        version     = "1.0",
        ext_modules = [_pyOrbfit]
        )

