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

where funcname is the name of the desired class or function.

"""

import numpy as _np
import vegas as _vegas
from .dukesMain import constant,snNuEenergy,_get_r,_dEv,vBDM,dmNumberDensity,FlagError,supernovaNuFlux
from .galDensity import galacticAreaDensity,galacticAreaDensityFit
from .galMassFunction import dnG,_E,rhoDotSFR


class userPhenoModelInterface(constant):
    """
    This class allows user to defined their own model-dependent differential cross sections of
    DM-nu and DM-e in CM frame.

    In principle, the differential cross sections depends on the Mandelstain variables (s,t,u)
    where they are composed of CM frame momentum, CM frame scattering angle and particle masses.
    The CM frame momentum can be expressed of incoming particle energy in the lab frame.


    //////////////////////////////////////
    //                                  //
    //          Initialization          //
    //                                  //
    //////////////////////////////////////
    
    To initialize the class, it takes two input functions func1 and func2 where func1 dictates
    the differential DM-nu cross section and func2 the differential DM-e cross section, both in
    CM frame. Note that each funciton can take only 3 arguments which explains below.

    ---*--- How to define your own cross sections ---*---
    
    For instance, one can define:

    # Differential DM-nu cross section
    def func1(Ev,mx,thetaCM_vx) -> float:
        ... body ...
        return some_value

    # Differential DM-e cross section
    def func2(Tx,mx,thetaCM_xe) -> float:
        ... body ...
        return some_value

    where some_val is the value of cross section in the unit of cm^2.

    You could have additional parameters in your pheno model such as coupling constant and
    mediator mass. These parameters are usually scalable and independent of the input arguments.
    So they are not necessary to be put in the function body unless they are running with the
    energy (running couplings) or mx-dependnet.


    ---*--- CAUTION ---*---
    
    If you designed func1 and 2 with more than 3 arguments, the surplus one(s) will be not be
    used in the flux and event calculations.


    ---*--- Initialized instance and check your implementation ---*---
    
    If one sucessfully initialized the userPhenoModel instance, the associated differential cross
    sections can be checked through

    # Initializing instance
    >>> myOwnModel = userPhenoModelInterface(func1,func2)
    
    # Input arguments
    >>> mx = 0.1
    >>> Ev,thetaCM_vx = 10,0.123
    >>> Tx,thetaCM_ve = 8.5,0.057
    
    # Differential cross sections of my model in the instance
    >>> myOwnModel.dsigmaNu(Ev,mx,thetaCM_vx)  # DM-nu, cm^2
    >>> myOwnModel.dsigmaE(Tx,mx,thetaCM_ve)   # DM-e, cm^2

    If they yield the values match your expectation, then the implementation is sucessful.


    ////////////////////////////////////////////
    //                                        //
    //   Running Flux and Event Calculation   //
    //                                        //
    ////////////////////////////////////////////

    For evaluating DBDM flux from model-dependent cross section, all inputs are the same as
    dukes.flux.

    For example,
    
    # Evaluating flux
    >>> Tx,mx = 5,1e-2
    >>> myOwnModel.flux(Tx=Tx,mx=mx)

    The output unit for flux is per MeV per cm^2 per second

    # Evaluating event
    >>> myOwnModel.event(mx=mx)

    The output unit for flux is per electron per second

    Though both func1 and func2 for model-dependent DM-nu and DM-e differential cross sections,
    respectively, are given during the initialization of the userPhenoModelInterface instance,
    func2 is not used in evaluating the flux as this does not require DM-e interaciton. func2
    is used when evaluating event.
    """

    def __init__(self,func1,func2):
        self.dsigmaNu = func1  # differential DM-nu cross section, takes 3 arguments: (Ev,mx,thetaCM_vx)
        self.dsigmaE = func2   # differential DM-e cross section, takes 3 arguments: (Tx,mx,thetaCM_xe)
    
    def _diffSpectrum(self,Tx,mx,MG,R,l,theta,thetaCM_vx,is_spike,sigv,tBH,rhosMW,rsMW,eta) -> float:
        """
        dNx/dTx
        """
        r = _get_r(l,R,theta)
        if r >= 1e-8:
            Ev = snNuEenergy(Tx,mx,thetaCM_vx)
            dEvdTx = _dEv(Tx,mx,thetaCM_vx)
            vx = vBDM(Tx,mx)  #  
            nx = dmNumberDensity(r,mx,MG,is_spike,sigv,tBH,rhosMW,rsMW,eta)
            dsigma0 = self.dsigmaNu(Ev,mx,thetaCM_vx) #self.dsigmaNu(Ev,mx,thetaCM_vx)
            return l**2*_np.sin(theta)*_np.sin(thetaCM_vx)*nx*dsigma0*supernovaNuFlux(Ev,l)*(dEvdTx*vx)
        else:
            return 0
    
    def _dbdmSpectrum(self,z,m,Tx,mx,R,l,theta,thetaCM_vx,is_spike,sigv,tBH,rhosMW,rsMW,eta) -> float:
        """
        DBDM spectrume yielded by SN at arbitrary position R
        """
        Txp = (1 + z)*Tx 
        if Txp < 200:  # discard the BDM signature if it requires Ev > 200 MeV at z 
            MG = 10**m
            return MG*dnG(m,z)/_E(z)*rhoDotSFR(z)*self._diffSpectrum(Txp,mx,MG,R,l,theta,thetaCM_vx,is_spike,sigv,
                                                                     tBH,rhosMW,rsMW,eta)
        else:
            return 0
        
    def _dbdmSpectrumWeighted(self,z,m,Tx,mx,R,l,theta,thetaCM_vx,is_spike,sigv,tBH,rhosMW,rsMW,eta,usefit) -> float:
        """
        DBDM spectrume yielded by SN at position R weighted by galactic baryonic distribution
        """
        Txp = (1 + z)*Tx 
        if Txp < 200:  # discard the BDM signature if it requires Ev > 200 MeV at z
            MG = 10**m
            # adopt fitting data for galactic area density?
            if usefit is True:
                galArealDensity = galacticAreaDensityFit((R,m))
            elif usefit is False:
                galArealDensity = galacticAreaDensity(R,zRange=[-10,10],MG=MG)
            else:
                raise FlagError('Global flag \'usefit\' must be a boolean.')
            
            return 2*_np.pi*R*galArealDensity*dnG(m,z)/_E(z)*rhoDotSFR(z)*self._diffSpectrum(Tx,mx,MG,R,l,theta,thetaCM_vx, 
                                                                                             is_spike,sigv,tBH,rhosMW,rsMW,eta)
        else:
            return 0

    def flux(self,Tx,mx,                                                           
             R=0,Rmax=500,rmax=500,tau=10,is_spike=True,is_average=True,      
             sigv=None,tBH=1e9,rhosMW=184,rsMW=24.42,eta=24.3856,usefit=True, 
             nitn=10,neval=50000) -> float:
        """
        DBDM flux for given (Tx,mx) for model-dependent DM-nu and DM-e differential cross
        sections
    
        In
        ------
        Tx: BDM kinetic energy, MeV
        mx: DM mass, MeV
        R: Specified SN position on the galactic plane, R=0=GC
            This only works for is_average = False
        Rmax: Maximum radius of galactic plane to be integrated
        rmax: Maximum halo radius to be integrated
        tau: SN duration, s
        is_spike: Including DM spike in the halo, bool
        is_average: SN position weighted by baryonic distribution, bool
        sigv: DM annihilation cross section, in the unit of 3e-26 cm^3/s
            None indicates no annihilation
        tBH: BH age
        rhosMW: NFW characteristic density for MW
        rsMW: NFW characteristic radius for MW
        eta: Mmw/Mhalo for MW
        nitn: Number of chains in vegas
        neval: Number of evaluation in each MCMC chain
    
        Out
        ------
        Flux: per MeV per cm^2 per second
        """
        preFactor = self.MagicalNumber  # constant.D_H0*0.017/constant.Mmw/rhoDotSFR(0)/1e6/constant.kpc2cm**2/constant.year2Seconds
        lmax = Rmax + rmax
        if is_average is True:
            integrator = _vegas.Integrator([[0,8],[6,12],[0,Rmax],[0,lmax],[0,_np.pi],[0,_np.pi]]) #(z,m,R,l,theta,thetaCM_vx)
            result = integrator(lambda x: self._dbdmSpectrumWeighted(z=x[0],m=x[1],Tx=Tx,mx=mx,R=x[2],l=x[3],theta=x[4],thetaCM_vx=x[5],    
                                                                     is_spike=is_spike,sigv=sigv,tBH=tBH,rhosMW=rhosMW,rsMW=rsMW,eta=eta,usefit=usefit),
                                nitn=nitn,neval=neval).mean
            flux = 4*_np.pi**2*tau*result*self.kpc2cm**3*vBDM(Tx,mx)*preFactor
        elif is_average is False:
            integrator = _vegas.Integrator([[0,8],[6,12],[0,lmax],[0,_np.pi],[0,_np.pi]]) #(z,m,l,theta,thetaCM_vx)
            result = integrator(lambda x: self._dbdmSpectrum(z=x[0],m=x[1],Tx=Tx,mx=mx,R=R,l=x[2],theta=x[3],thetaCM_vx=x[4], 
                                                             is_spike=is_spike,sigv=sigv,tBH=tBH,rhosMW=rhosMW,rsMW=rsMW,eta=eta),      
                                nitn=nitn,neval=neval).mean
            flux = 4*_np.pi**2*tau*result*self.kpc2cm**3*vBDM(Tx,mx)*preFactor
        else:
            raise FlagError('Flag \'is_average\' must be a boolean.')
        return flux

    def event(self,mx,                                                                                                                
          TxRange=[5,100],R=0,Rmax=500,rmax=500,tau=10,is_spike=True,is_average=True,  
          sigv=None,tBH=1e9,rhosMW=184,rsMW=24.42,eta=24.3856,usefit=True,             
          nitn=10,neval=50000) -> float:
        """
        DBDM event per electron per second for given mx for model-dependent DM-nu and DM-e
        differential cross sections 
    
        In
        ------
        mx: DM mass, MeV
        TxRange: [Tx_min,Tx_max], Tx range to be integrated over, MeV
        R: Specified SN position on the galactic plane, R=0=GC
            This only works for is_average = False
        Rmax: Maximum radius of galactic plane to be integrated
        rmax: Maximum halo radius to be integrated
        tau: SN duration, s
        is_spike: Including DM spike in the halo, bool
        is_average: SN position weighted by baryonic distribution, bool
        sigv: DM annihilation cross section, in the unit of 3e-26 cm^3/s
            None indicates no annihilation
        tBH: BH age
        rhosMW: NFW characteristic density for MW
        rsMW: NFW characteristic radius for MW
        eta: Mmw/Mhalo for MW
        nitn: Number of chains in vegas
        neval: Number of evaluation in each MCMC chain
    
        Out
        ------
        Event: per electron per second
        """
        preFactor = self.MagicalNumber  # constant.D_H0*0.017/constant.Mmw/rhoDotSFR(0)/1e6/constant.kpc2cm**2/constant.year2Seconds
        lmax = Rmax + rmax
        if is_average is True:
            integrator = _vegas.Integrator([[0,8],[6,12],[0,Rmax],[0,lmax],[0,_np.pi],[0,_np.pi],[0,_np.pi],TxRange]) #(z,m,R,l,theta,thetaCM_vx,thetaCM_xe,Tx)
            result = integrator(lambda x: self._dbdmSpectrumWeighted(z=x[0],m=x[1],Tx=x[7],mx=mx,R=x[2],l=x[3],theta=x[4],thetaCM_vx=x[5],  
                                               is_spike=is_spike,sigv=sigv,tBH=tBH,         
                                               rhosMW=rhosMW,rsMW=rsMW,eta=eta,usefit=usefit)*vBDM(Tx=x[7],mx=mx)*self.dsigmaE(x[7],mx,x[6])*_np.sin(x[6]), 
                            nitn=nitn,neval=neval).mean
            event = 4*_np.pi**2*tau*result*self.kpc2cm**3*preFactor*(2*_np.pi)  # the last 2pi is due to diff cross section for DM-e is independent of azimuthal angle 
        elif is_average is False:
            integrator = _vegas.Integrator([[0,8],[6,12],[0,lmax],[0,_np.pi],[0,_np.pi],[0,_np.pi],TxRange]) #(z,m,l,theta,thetaCM,thetaCM_xe,Tx)
            result = integrator(lambda x: self._dbdmSpectrum(z=x[0],m=x[1],Tx=x[6],mx=mx,R=R,l=x[2],theta=x[3],thetaCM_vx=x[4],     
                                               is_spike=is_spike,sigv=sigv,tBH=tBH,         
                                               rhosMW=rhosMW,rsMW=rsMW,eta=eta,usefit=usefit)*vBDM(Tx=x[6],mx=mx)*self.dsigmaE(x[6],mx,x[5])*_np.sin(x[5]), 
                                nitn=nitn,neval=neval).mean
            event = 4*_np.pi**2*tau*result*self.kpc2cm**3*preFactor*(2*_np.pi)  # the last 2pi is due to diff cross section for DM-e is independent of azimuthal angle 
        else:
            raise FlagError('Flag \'is_average\' must be a boolean.')
        return event
