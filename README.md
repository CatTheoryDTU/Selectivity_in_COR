# Data and analysis tools for the study of product selectivity in eCO2R
# published in ACS Catal. 2023, 13, 7, 5062–5072 (https://pubs.acs.org/doi/full/10.1021/acscatal.3c00228)

The repository contains the bare DFT optimized geometries and energies included
in the folder `dft_data`. Here the thermodynamics (stable intermediates) and
kinetc minimum energy paths are given in separate directories. All the data has
been calculated on Cu(100).

Furthermore, the routine to reproduce the figures in the paper are given in the
directory `figures`. Each directory contains its own README.

All the calculations and analysis have been performed using the atomic
simulation environment (ASE) version 3.19.0b1 and GPAW version 21.1.1b1. The
microkinetic model was simulated using CatMap from
https://github.com/sringe/catmap-1.

WARNING: Sadly, in order to parse the trajectory files properly a slightly changed ASE version from trunk is needed, as vanilla ASE does not allow readinig and writing of the number of electrons and potentials in the trajectory format. The necessary changes to ASE are:

In ase/calculators/singlepoint.py change line 25 from
```
if property in ['energy', 'magmom', 'free_energy']:
```
to
```
if property in ['energy', 'magmom', 'free_energy','ne','electrode_potential']:
```


and in ase/calculators/calculator line 98 exchange:
```
all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy']
```
with
```
all_properties = ['energy', 'forces', 'stress', 'stresses', 'dipole',
                  'charges', 'magmom', 'magmoms', 'free_energy','ne','electrode_potential']
```

We also note that newer GPAW versions changed the names of `ne` and `electrode_potential`. However, ASE can still read the results if above's changes are made.


