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

__name__         = 'dukes'
__version__      = '1.2.0'
__description__  = 'This package evaluates the signatures of diffuse boosted dark matter by supernova neutrinos in the early Universe'
__author__       = 'Yen-Hsun Lin'
__email__        = 'yenhsun@phys.ncku.edu.tw'
__url__          = 'https://github.com/yenhsunlin/dukes'
__license__      = 'GNU GPL-3.0'

from .dukesMain import *
from .galDensity import generalDensityProfile,mwDensityProfile,galacticDensityProfile
