# `doom`: **D**iffuse b**OO**sted dark **M**atter yielded from supernova neutrinos in the early universe


## Introduction

`doom` is a package for evaluating the diffuse signature of dark matter boosted in the early Universe due to supernova neutrinos based on `arXiv:20xx.xxxx`.

### Installation

To install, excute following command on the terminal:

    $ pip install doom

and everything should be processed on-the-fly.

### Dependencies

`doom` requires these external packages

- `numpy` >= 1.20.0
- `scipy` >= 1.10.0
- `vegas` >= 6.0.1

where `vegas` is a the backend engine for evaluating multidimensional integrals based on adaptive Monte Carlo vegas algorithm, see its homepage: [https://pypi.org/project/vegas/](https://pypi.org/project/vegas/).

Other packages, e.g. `gvar`, maybe required by these dependencies during the installation.
The versions of these dependencies are not strict, but are recommended to update to the latest ones to avoid incompatibility. 


## Usage

We briefly summarize the usage in this section and a comprehensive tutorial can be found in the jupyter notebook stored in `tutorial/tutorial.ipynb`.

To import, do

    >>> import doom

in the python terminal and is similar in the jupyter notebook. All module functions can be called like `dbdm.funcname`.

### Examples

#### Physical constants

#### Boosted dark matter velocity

A boosted dark matter (BDM) with mass $m_\chi$ and kinetic energy $T_\chi$ has the velocity $v_\chi$,

$$
\frac{v_\chi}{c} = \frac{\sqrt{T_\chi(2m_\chi+T_\chi)}}{m_\chi+T_\chi}.
$$

Let $T_\chi=$ `Tx` and $=m_\chi=$ `mx`, the corresponding function that evaluates $v_\chi/c$ is `doom.vBDM(Tx,mx)` is

    >>> Tx,mx = 5,1  # MeV
    >>> doom.vBDM(Tx,mx)
    0.9860132971832694
## Misc

Bug report and troubleshooting please contact the author Yen-Hsun Lin via [yenhsun@phys.ncku.edu.tw](mailto:yenhsun@phys.ncku.edu.tw).
