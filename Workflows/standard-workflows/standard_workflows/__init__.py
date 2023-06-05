#!/usr/bin/env python

#
# Copyright 2022, Saez Lab
#
# File author(s): Hanna Schumacher <hanna@schumacher-home.de>
#                 Denes Turei <turei.denes@gmail.com>
#
# Distributed under GPLv3 license, see the file `LICENSE`.
#


"""
Single Cell Standard Analysis Framework
"""

#__all__ = ['twentythree']

from ._metadata import __author__, __version__
from . import memoize
from . import sc_analysis_baseclass as sc_classes
from . import sc_analysis_loops as scl
from . import sc_decoupler_utility as dcu
from . import scfunctions as sc_funcs
import yaml, json

def twentythree() -> int:
    """
    The number twenty-three.
    """
   
    return 23

