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


class constant:
    """
    Containing multiple physical constants and conversion factors
    
    Conversion factors
    ------
    md2MeVperCubicCM: Msun/kpc^3 to MeV/cm^3
    md2MeVperQuadCM: Msun/kpc^3 to MeV/cm^2
    year2Seconds
    erg2MeV
    
    Masses (MeV)
    ------
    me: electron
    mn: neutron
    mp: proton
    Msun: Solar mass, MeV
    Msun_kg: Solar mass, kg
    Mmw: Milky Way stellar mass, Msun
    Mhalo: DM halo mass of MW, Msun
    
    Lengths (cm)
    ------
    Rhalo: estimated MW halo radius, kpc

    Cross section
    ------
    sigma0: default value for cross section, cm^2
    
    Physical constants
    ------
    c: light speed, cm/s
    H0: Hubble constant, km/s/Mpc
    rho_c: critical density, Msun/pc^3
    Lv: SN neutrino luminosity, single specie, erg/s
    G: Newton's constant, m^3/kg/s^2
    """
    # conversion factors 
    md2MeVperCubicCM = 3.796e-05        # Msun/kpc^3 to MeV/cm^3
    md2MeVperQuadCM  = 1.171e+17        # Msun/kpc^3 to MeV/cm^2
    year2Seconds     = 3.156e+07        # year to seconds
    erg2MeV          = 6.241e+05        # erg to MeV
    kpc2cm           = 3.085e+21        # kpc to cm
    
    # masses
    me               = 5.110e-01        # electron mass, MeV
    mn               = 9.395e+02        # neutron mass, MeV
    mp               = 9.382e+02        # proton mass, MeV
    Msun             = 1.115e+60        # Solar mass, MeV
    Msun_kg          = 1.981e+30        # Solar mass, kg
    Mmw              = 5.290e+10        # MW stellar mass, Msun
    Mhalo            = 1.290e+12        # MW halo mass, Msun
    
    # lengths
    Rhalo            = 2.300e+02        # MW halo radius, kpc

    # cross section
    sigma0           = 1.000e-35        # sigma_0, cm^2

    # physical constants
    c                = 2.998e+10        # light speed, cm/s
    H0               = 7.300e+10        # Hubble constant, km/s/Mpc
    rho_c            = 1.500e-07        # critical density, Msun/pc^3
    Lv               = 3.000e+52/6      # Supernova neutrino luminosity (single specie), erg/s
    Omega_0m         = 3.150e-01        # Cosmological matter fraction
    Omega_0L         = 6.850e-01        # Cosmological dark energy fraction
    Omega_0r         = 2.300e-03        # Cosmological radiation fraction
    Omega_0          = 1.000            # Cosmological total energy
    D_H0             = 4.280e+03        # Mpc
    G                = 6.674e-11        # Newton gravitational constant, meter^3 kg^-1 s^-2

    MagicalNumber    = 2.572e-64