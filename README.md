[![Python](https://img.shields.io/badge/python-3.8-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/License-GPL_3.0-blue.svg)](https://choosealicense.com/licenses/gpl-3.0/)
[![ArXiv](https://img.shields.io/badge/arXiv-2404.08528-yellowgreen.svg)](https://arxiv.org/abs/2404.08528) 


# dukes: *D*iff*U*se-boosted dar*K* matt*E*r by *S*upernova neutrinos


`dukes` is a pyhon package for evaluating the signatures of diffuse boosted dark matter by supernova neutrinos in the early Universe based on [arXiv:2404.08528](https://arxiv.org/abs/2404.08528).
It also supports an *experimental* feature that implements particle-model-dependent differential cross sections for DM-neutrino and DM-electron.


### Citation

If you use this package or part of the code in your research, please cite the following:

1. Yen-Hsun Lin and Meng-Ru Wu, *Echoes of darkness: Supernova-neutrino-boosted dark matter from all galaxies*, arXiv:2404.08528
2. `dukes`: https://github.com/yenhsunlin/dukes

## <span style="color:cyan">SCHEDULED MAJOR UPDATE OF VERSION 2</span>

Through the release of `snorer`, a python package for evaluating time-of-flight SN*ν* BDM signature in Milky Way, Large Magellanic Cloud and arbitrary distant galaxy, the functions and classes for dealing with DM halo shape, BDM kinematics and various useful astrophysical relations become mature, clear and reusable. Furthermore, the halo spike for $\alpha=7/3$ is also included with custmizable halo shape for users.

In the next major update of `dukes`, we will fully import these new features from `snorer` and the old ones will be replaced or depricated.
The two versions are generall *incompatible*.
The changes and their new features are listed below:

|  original   | replaced/added  |type|new features|
|  ----  | ----  | ---- |---- |
| `haloSpike`  | `HaloSpike` |*class*|including $\alpha=7/3$ and customizable halo shape|
||`Kinematics`|*class*|Handling  kinematics of $2\to2$ particle scattering|
||`Mandelstam`|*class*|Handling Mandelstam variables|
| `constant` | `Constants` |*class*|Various new data added|
||`constant`|*inst*|Instance of `Constants`|
|`FlagError`|`FlagError`|*class*|Repeated feature|
|`dmNumberDensity()`|`dmNumberDensity()`|*func*|Customizable halo shape|
|`snNuEenergy()`|`snNuSpectrum()`|*func*|Allowing number density output and truncation point added|
|`rhox()`|`rhox()`|*func*|Allowing customizable halo slope $n$|
|`nxNFW()`||*func*|Fully integrated into `dmNumberDensity()`|
|`radiusSchwarzschild()`|`radiusSchwarzschild()`|*func*|Repeated feature|
||`M_sigma()`|*func*|$M-\sigma$ relation|
||`radiusInfluence()`|*func*|SMBH influence radius|

If there is no original one, it means that one is newly added.
On the other hand, if there is no replaced one, it implies that one will be depricated.
Even for functions and classes that are not being replaced, their features might be upgraded as well. We will indicate them in the new `tutorial.ipynb`.

## Installation

To install, excute the following command on the prompt

    $ pip install dukes

and everything should be processed on-the-fly.

### Dependency

`dukes` requires python >= 3.8 and the following packages

- `numpy` >= 1.20.0
- `scipy` >= 1.10.0
- `vegas` >= 6.0.1

where `vegas` is the backend engine for evaluating multidimensional integrals based on adaptive Monte Carlo vegas algorithm, see its homepage: [https://pypi.org/project/vegas/](https://pypi.org/project/vegas/).

Other packages, e.g. `gvar`, maybe required by these dependencies during the installation.
The versions of these dependencies are not strict, but are recommended to update to the latest ones to avoid incompatibility. 


## Usage

We briefly summarize the usage in this section and a comprehensive tutorial can be found in the jupyter notebook in `examples/tutorial.ipynb`.

To import, do

    >>> import dukes

in python terminal and is similar in the jupyter notebook. All module functions named *funcname* can be called by typing `dukes.funcname`.


### Boosted dark matter velocity

A boosted dark matter (BDM) with mass $m_\chi$ and kinetic energy $T_\chi$ has the velocity $v_\chi$,

$$
\frac{v_\chi}{c} = \frac{\sqrt{T_\chi(2m_\chi+T_\chi)}}{m_\chi+T_\chi}.
$$

Let $T_\chi=$ `Tx` and $=m_\chi=$ `mx`, the corresponding function that evaluates $v_\chi/c$ is

    >>> Tx,mx = 5,1  # MeV
    >>> dukes.vBDM(Tx,mx)
    0.9860132971832694


### The diffuse BDM flux

The averaged diffuse BDM (DBDM) flux at redshift $z=0$ is given by

$$
\frac{d\Phi_\chi}{dT_\chi} = \frac{v_\chi}{H_0} \int_0^{z_{\rm max}} \frac{dz}{\varepsilon(z)}  \int dM_G \frac{d\Gamma_{{\rm SN}}(z)}{dM_G}\frac{d\bar N_\chi(M_G)}{dT_\chi^\prime}. 
$$

Same as the above example,

    >>> dukes.flux(Tx,mx,usefit=True,nitn=10,neval=50000)
    4.422705310516041e-08

as in MeV<sup>−1</sup> cm<sup>−2</sup> s<sup>−1</sup>.

Throughout the entire package, we have implemented the differential DM-neutrino scattering cross section in CM frame is isotropic and energy-independent

$$
\frac{d\sigma_{\chi \nu}}{d\Omega_{\rm CM}}=10^{-35}~{\rm cm^2~sr^{-1}}.
$$

The argument `usefit` is to turn on/off the fitting function used in obtaining the average supernova position on the galactic plane for galaxy with baryonic mass $M_G$.
If `usefit=False`, the function will call `galacticDensityProfile` to evaluate the area density for galaxy with arbitrary $M_G$.
It requires quadrature integration `quad` from scipy and the computation time surges accordingly, but the accuracy is improved insignificantly.

The arguments `nitn` and `neval` are passed to `vegas` and determine how many chains of iteration and how many numbers to be evaluated in each chain. Increasing them will improve the accuracy of the results but also cost longer computation time. We relegate the detail to `vegas` [documentation](https://vegas.readthedocs.io/).


### Physical constants

We have a class named `constant` that contains multiple physical constants and conversion factors frequently used in this package.
For instance, electron mass

    >>> dukes.constant.me
    0.511

as in MeV and the speed of light

    >>> dukes.constant.c
    29980000000.0

as in cm s<sup>−1</sup>.
Conversion factors such as converting kiloparsec to centimeters

    >>> dukes.constant.kpc2cm
    3.085e21

and year to seconds

    >>> dukes.constant.year2Seconds
    31556926


## Scripting

In python script (see `tests/dukes_example.py`), one can write

    # dukes_example.py

    import sys
    import dukes

    if __name__ == '__main__':

        Tx = float(sys.argv[1])  # DM kinetic energy, MeV
        mx = float(sys.argv[2])  # DM mass, MeV
        vx = dukes.vBDM(Tx,mx)   # BDM velocity
        
        print(vx)                # Print the BDM velocity

and excute this on the prompt

    $ python dukes_example.py 5 1
    0.9860132971832694

or whatever style you like!

## Bugs and troubleshooting

Please report to the author, Yen-Hsun Lin, via [yenhsun@phys.ncku.edu.tw](mailto:yenhsun@phys.ncku.edu.tw).