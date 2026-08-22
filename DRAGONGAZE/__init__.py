'''
DRAGONGAZE
Created by Bobbie Markwick.

Marks this directory as a Python package so its modules can use package-relative
imports (`from .Config import *`, etc). This is what lets Config.py/HDF5_keys.py/
FileUtils.py be shared unambiguously between the GAZE*.py scripts without relying
on the current working directory or sys.path tricks.

Consequence: because of the relative imports, the GAZE*.py scripts can no longer
be run directly (`python GAZE1D.py`) from inside DRAGONGAZE/ -- Python raises
"attempted relative import with no known parent package" for a script executed
that way, since it isn't run as part of a package. Instead, run them as modules
from the DIRECTORY THAT CONTAINS DRAGONGAZE/ (e.g. .../Code/DRAGON/), e.g.:

    python -m DRAGONGAZE.GAZE2D rho p Bmag
'''
