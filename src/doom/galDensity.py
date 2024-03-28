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


#from os.path import dirname
import numpy as _np
from scipy.interpolate import CubicSpline as _CubicSpline
from scipy.interpolate import RegularGridInterpolator as _RegularGridInterpolator
from scipy.integrate import quad as _quad
from .dat.densityParamFit import _MG_density_data,_rho_b0_data,_Sigma0Thick_data,_Sigma0Thin_data
from .dat.galcticAreaDensityFit import _R_data,_MG_area_data,_rho_data



##########################################################################
#                                                                        #
#   Fitting Functions from Data Points                                   #
#                                                                        #
##########################################################################


# ----- Load data points ----- 
# Data points for parameters in describing baryonic distribution for arbitrary
# galactic mass MG. These parameters are scaled from Milky Way case documented
# in MNRAS 465, 76 (2017).

#_densityParamData = _np.load(dirname(__file__) + '/bin/densityParamFit.npz')
#_rho_b0_data = _densityParamData['rho_b0']            # rho_b0 for bulge, log_10
#_Sigma0Thin_data = _densityParamData['Sigma0Thin']    # Sigma_0 for stellar thin, log_10
#_Sigma0Thick_data = _densityParamData['Sigma0Thick']  # Sigma_0 for stellar thick, log_10 
#_MG_density_data = _densityParamData['MG']            # Galactic mass MG, log_10

# Data points for baryonic area density at given location R for arbitrary
# galactic mass MG.
#_areaDensityData = _np.load(dirname(__file__) + '/bin/galcticAreaDensityFit.npz')
#_R_data = _areaDensityData['R']                       # Position, kpc
#_MG_area_data = _areaDensityData['MG']                # Galactic mass MG, log_10
#_rho_data = _areaDensityData['AreaDensity']           # Area density at (R,MG)


# ----- Fitting data points ----- 
_rho_b0Fit = _CubicSpline(_MG_density_data,_rho_b0_data)
_Sigma0ThinFit = _CubicSpline(_MG_density_data,_Sigma0Thin_data)
_Sigma0ThickFit = _CubicSpline(_MG_density_data,_Sigma0Thick_data)
galacticAreaDensityFit = _RegularGridInterpolator((_R_data,_MG_area_data),_rho_data)



##########################################################################
#                                                                        #
#   General Classes and Functions for Numerics                           #
#                                                                        #
##########################################################################


class generalDensityProfile:
    """
    This class contains the general mathematical expressions for bulge, stellar
    disc and gas disc based on McMillan, MNRAS 465, 76 (2017). The unit is M_sun
    per cubic kpc.
    """
    
    def __init__(self):
        pass
    
    @staticmethod
    def bulgeForm(R,z,rho_0b,r_0,r_cut,a,q) -> float:
        """
        Eq. (1)
        """
        r_p = _np.sqrt(R**2 + (z/q)**2)
        rho_b = rho_0b*_np.exp(-r_p/r_cut)/(1 + r_p/r_0)**a
        return rho_b
    
    @staticmethod
    def stellarForm(R,z,Rd,zd,Sigma0) -> float:
        """
        Eq. (3)
        """
        return Sigma0*_np.exp(-_np.abs(z)/zd - R/Rd)/2/zd
    
    @staticmethod
    def gasForm(R,z,Rd,Rm,zd,Sigma0) -> float:
        """
        Eq. (4)
        """
        R = _np.abs(R)
        return Sigma0*_np.exp(-Rm/R - R/Rd)/_np.cosh(z/2/zd)**2/4/zd


def mwDensityProfile(R,z) -> float:
    """
    Milky Way density profile at given position
    
    In
    ------
    R: Distance to the MW center, kpc
    z: Height to the galactic plane, kpc
        Negative value means below the plane
    
    Out
    ------
    density: M_sun/kpc^3
    """
    density = 0
    density += generalDensityProfile.stellarForm(R,z,2.5,0.3,8.96e8)           # stellar thin
    density += generalDensityProfile.stellarForm(R,z,3.02,0.9,1.83e8)          # stellar thick
    #density += generalDensityProfile.gasForm(R,z,7,4,0.085,5.31e7)             # gas HI
    #density += generalDensityProfile.gasForm(R,z,1.5,12,0.045,2.18e9)          # gas H2
    density += generalDensityProfile.bulgeForm(R,z,9.93e10,0.075,2.1,1.8,0.5)  # bulge
    return density


def galacticDensityProfile(R,z,MG) -> float:
    """
    Generalized density profile with arbitrary galactic mass MG at given position
    
    In
    ------
    R: Distance to the galactic center, kpc
    z: Height to the galactic plane, kpc
        Negative value means below the plane
    MG: The galactic mass, Msun
    
    Out
    ------
    density: M_sun/kpc^3
    """
    # Here we assume that the ratios of bulge mass and stellar mass to MG are
    # the same as Milky Way for all galaxies for simplicity. In addition, we
    # also assume the stellar thin and thick ratio is the same too
    # Mb:Ms = 0.163:0.837 and M_thin:M_thick = 0.77:0.23
    massBulge = 0.163*MG   
    massStellarThin = 0.837*0.77*MG
    massStellarThick = 0.837*0.23*MG 
    
    # parameters for bulge
    a = 1.8
    q = 0.5
    r_0 = 0.075
    rho_0b = 10**_rho_b0Fit(_np.log10(massBulge))
    r_cut = 2.1*(massBulge/8.958e9)**(1/3)
    r_p = _np.sqrt(R**2 + (z/q)**2)
    rho_b = rho_0b*_np.exp(-r_p/r_cut)/(1 + r_p/r_0)**a
    
    # parameters for stellar thin
    zd_thin = 0.3*(massStellarThin/3.518e10)**0.5
    Rd_thin = 2.5*(massStellarThin/3.518e10)**0.5
    Sigma0_thin = 10**_Sigma0ThinFit(_np.log10(massStellarThin))
    
    # parameters for stellar thick
    zd_thick = 0.9*(massStellarThick/1.048e10)**0.5
    Rd_thick = 3.02*(massStellarThick/1.048e10)**0.5
    Sigma0_thick = 10**_Sigma0ThickFit(_np.log10(massStellarThick))
    
    # calculate density
    density = 0
    density += generalDensityProfile.bulgeForm(R,z,rho_0b,r_0,r_cut,a,q)
    density += generalDensityProfile.stellarForm(R,z,Rd_thin,zd_thin,Sigma0_thin)
    density += generalDensityProfile.stellarForm(R,z,Rd_thick,zd_thick,Sigma0_thick)
    
    return density


def galacticAreaDensity(R,zRange=[-10,10],MG=None) -> float:
    """
    Evaluates the area density at given R after height is integrated out.
    
    In
    ------
    R: Distance to the galactic center, kpc
    zRange: The range of z to be integrated out, Default
        is [-10,10]
    MG: The galactic mass, Msun
    
    Out
    ------
    area density: M_sun/kpc^2
    """
    zmin,zmax = zRange
    if MG is None:
        return _quad(lambda z: mwDensityProfile(R,z),zmin,zmax)[0]
    else:
        return _quad(lambda z: galacticDensityProfile(R,z,MG),zmin,zmax)[0]