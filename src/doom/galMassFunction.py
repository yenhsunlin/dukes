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


import numpy as _np
from .constant import constant



##########################################################################
#                                                                        #
#   General Functions for Numerics                                       #
#                                                                        #
##########################################################################


def dnG(m,z):
    """
    Stellar mass function, dnG/dm.
    
    In
    ------
    m: Log10(MG/Msun)
    z: Redshift
    
    Out
    ------
    # of galaxies per m per Mpc^3
    """
    
    def phi(alpha,Mc,phi0,M):
        log_phi0 = _np.log10(phi0) - 4
        logln10 = 0.362216
        logE = 0.43429
        logphi = log_phi0 + logln10 + (M - Mc)*(1 + alpha) - 10**(M - Mc)*logE
        return 10**logphi
    
    if 0 <= z <= 0.7:
        return phi(-1.11,11.22,18.2,m)
    elif 0.7 < z <= 1:
        return phi(-1.27,11.37,11.0,m)
    elif 1 < z <= 1.4:
        return phi(-1.28,11.26,6.2,m)
    elif 1.4 < z <= 1.8:
        return phi(-1.31,11.25,4.3,m)
    elif 1.8 < z <= 2.2:
        return phi(-1.34,11.22,3.1,m)
    elif 2.2 < z <= 2.6:
        return phi(-1.38,11.16,2.4,m)
    elif 2.6 < z <= 3:
        return phi(-1.41,11.09,1.9,m)
    elif 3 < z <= 3.5:
        return phi(-1.45,10.97,1.5,m)
    elif 3.5 < z <= 4:
        return phi(-1.49,10.81,1.1,m)
    elif 4 < z <= 4.5:
        return phi(-1.53,10.44,3.0,m)
    elif 4.5 < z <= 5.5:
        return phi(-1.67,10.47,1.3,m)
    elif 5.5 < z <= 6.5:
        return phi(-1.93,10.3,0.3,m)
    elif 6.5 < z <= 8:
        return phi(-2.05,10.42,0.1,m)
    else:
        raise ValueError('Argument \'z\' is not in the valid range.')
    

def _E(z):
    y = (1 + z)
    return _np.sqrt(constant.Omega_0m*y**3 + constant.Omega_0L)


def rhoDotSFR(z):
    """
    Star formation rate at redshift z
    
    In
    ------
    z: Redshift
    
    Out
    ------
    SFR rate: Msun per year per Mpc^3
    """
    a,b,c = 3.4,-0.3,-3.5
    z1,z2,eta = 1,4,-10
    B = (1 + z1)**(1 - a/b)
    C = (1 + z1)**((b - a)/c)*(1 + z2)**(1 - b/c)
    Z = 1 + z
    return 0.0178*(Z**(a*eta) + (Z/B)**(b*eta) + (Z/C)**(c*eta))**(1/eta)
