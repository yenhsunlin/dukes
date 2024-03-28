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
    Msun: Solar mass
    
    Lengths (cm)
    ------
    kpc2cm
    
    Physical constants
    ------
    c: light speed, cm/s
    H0: Hubble constant, km/s/Mpc
    rho_c: critical density, Msun/pc^3
    Lv: SN neutrino luminosity, single specie, erg/s
    """
    # conversion factors 
    md2MeVperCubicCM = 3.796e-5  # Msun/kpc^3 to MeV/cm^3
    md2MeVperQuadCM = 1.171e17   # Msun/kpc^3 to MeV/cm^2
    year2Seconds = 3.156e7       # year to seconds
    erg2MeV = 624150.913         # erg to MeV
    kpc2cm = 3.085e21            # kpc to cm
    
    # masses
    me = 0.511                   # electron mass, MeV
    mn = 939.564                 # neutron mass, MeV
    mp = 938.272                 # proton mass, MeV
    Msun = 1.115e60              # Solar mass, MeV
    Msun_kg = 1.981e30           # Solar mass, kg
    Mmw = 5.29e10                # MW stellar mass, Msun
    Mhalo = 1.29e12              # MW halo mass, Msun
    
    # lengths
    Rhalo = 230                  # MW halo radius, kpc

    # physical constants
    c = 2.998e10                 # light speed, cm/s
    H0 = 73                      # Hubble constant, km/s/Mpc
    rho_c = 1.5e-7               # critical density, Msun/pc^3
    Lv = 3e52/6                  # Supernova neutrino luminosity (single specie), erg/s
    Omega_0m = 0.315             # Cosmological matter fraction
    Omega_0L = 0.685             # Cosmological dark energy fraction
    Omega_0r = 0.0023            # Cosmological radiation fraction
    Omega_0 = 1                  # Cosmological total energy
    D_H0 = 4280                  # Mpc
    G = 6.6743e-11               # Newton gravitational constant, meter^3 kg^-1 s^-2

    MagicalNumber = 2.57258e-64