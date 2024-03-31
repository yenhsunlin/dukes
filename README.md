# dukes: *D*iff*U*se-boosted dar*K* matt*E*r by *S*upernova neutrinos


`dukes` is a package for evaluating the signatures of diffuse boosted dark matter by supernova neutrinos in the early Universe based on arXiv:24xx.xxxxx.
`dukes` also supports implementation of phenomenlogical-model-dependent differential cross sections between DM-neutrino and DM-electron for calculating BDM signatures (*experimental*).

## Installation

To install, excute the following command on the prompt

    $ pip install dukes

and everything should be processed on-the-fly.

### Dependency

`dukes` requires these external packages

- `numpy` >= 1.20.0
- `scipy` >= 1.10.0
- `vegas` >= 6.0.1

where `vegas` is a the backend engine for evaluating multidimensional integrals based on adaptive Monte Carlo vegas algorithm, see its homepage: [https://pypi.org/project/vegas/](https://pypi.org/project/vegas/).

Other packages, e.g. `gvar`, maybe required by these dependencies during the installation.
The versions of these dependencies are not strict, but are recommended to update to the latest ones to avoid incompatibility. 


## Usage

We briefly summarize the usage in this section and a comprehensive tutorial can be found in the jupyter notebook stored in `tutorial/tutorial.ipynb`.

To import, do

    >>> import dukes

in python terminal and is similar in the jupyter notebook. All module functions named *funcname* can be called by typing `dukes.funcname`.

### Examples

#### Boosted dark matter velocity

A boosted dark matter (BDM) with mass $m_\chi$ and kinetic energy $T_\chi$ has the velocity $v_\chi$,

$$
\frac{v_\chi}{c} = \frac{\sqrt{T_\chi(2m_\chi+T_\chi)}}{m_\chi+T_\chi}.
$$

Let $T_\chi=$ `Tx` and $=m_\chi=$ `mx`, the corresponding function that evaluates $v_\chi/c$ is

    >>> Tx,mx = 5,1  # MeV
    >>> dukes.vBDM(Tx,mx)
    0.9860132971832694


#### The diffuse BDM flux

The averaged diffuse BDM (DBDM) flux on the Earth is given by

$$
\frac{d\Phi_\chi}{dT_\chi} = \frac{v_\chi}{H_0} \int_0^{z_{\rm max}} \frac{dz}{\varepsilon(z)}  \int dM_G \frac{d\Gamma_{{\rm SN}}(z)}{dM_G}\frac{d\bar N_\chi(M_G)}{dT_\chi^\prime}. 
$$

Same as the above example,

    >>> dukes.flux(Tx,mx,usefit=True,nitn=10,neval=50000)
    4.422705310516041e-08

as in MeV<sup>−1</sup> cm<sup>−2</sup> s<sup>−1</sup>.

Throughout the entire package, we have implemented the differential DM-neutrino scattering cross section in CM frame is isotropic and energy-independent

$$
\frac{d\sigma_{\chi \nu}}{d\Omega_{\rm CM}}=\frac{\sigma_0}{4\pi}
$$

where $\sigma_0=10^{-35}$ cm<sup>2</sup>.

The argument `usefit` is to turn on/off the fitting function used in obtaining the average supernova position on the galactic plane for galaxy with baryonic mass $M_G$.
If `usefit=False`, the function will call `galacticDensityProfile` to evaluate the area density for galaxy with arbitrary $M_G$.
It requires quadrature integration `quad` from scipy and the computation time surges accordingly, but the accuracy is improved insignificantly.

The arguments `nitn` and `neval` are passed to `vegas` and determine how many chains of iteration and how many numbers to be evaluated in each chain. Increasing them will improve the accuracy of the results but also cost longer computation time. We relegate the detail to `vegas` [documentation](https://vegas.readthedocs.io/).


#### Physical constants

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
    31560000.0


### Scripting

In python script (see subsidiary `tests/dukes_example.py`), one can write

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