# Created by Yen-Hsun Lin (Academia Sinica) in 03/2024.
# Copyright (c) 2024 Yen-Hsun Lin.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version (see <http://www.gnu.org/licenses/>).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.



"""

This module contains utilities that are not necessary for the main
program. All the classes and functions in this module can be loaded
manually by users via

>>> from dukes.utils import funcname

where funcname are the desired class or function.

"""

import numpy as _np
import vegas as _vegas
from .dukesMain import constant,snNuEenergy,_get_r,_dEv,vBDM,dmNumberDensity
from .galDensity import galacticAreaDensity,galacticAreaDensityFit
from .galMassFunction import dnG,_E,rhoDotSFR

