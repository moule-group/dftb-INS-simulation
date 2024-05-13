# DFTB-INS-simulation

Use Density functional tight binding theory to simulate inelastic neutron scattering spectra.

## Installing DFTB 
Create the environment in your PC.
```
conda create -n dftb+
```
Or create the environment in supercomputer (ex. NERSC).
```
conda create -p /global/common/software/xxxx/dftb+
```
Install DFTB by using mamba.
```
conda activate dftb+
mamba install 'dftbplus=*=nompi_*'
```
Download Slater-Koster files from DFTB website.
```
https://dftb.org/parameters/download/all-sk-files
```
If you want to access DFTB by ASE, you have to set environment variable in your .bashrc file.
```
export DFTB_PREFIX="/path_to/slako/mio/mio-1-1/"
export ASE_DFTB_COMMAND="/path_to/dftb+/bin/dftb+ > PREFIX.out"
```
Add another environment variable in your .bashrc file in order to run dftb+ without activating envs everytime.
```
export PATH"=/path_to/dftb+/bin:$PATH"
```
Now, we finish the installation of DFTB! :smirk:

## Simulated INS spectras by DFTB (dftb.py)
Main functions:
```
runDFTB(kPoints,geometry=None,temp=5,fmax=0.01,mode=None)
```
It will automatically find the structure file in the folder(*.cif, *.gen)

There are two modes, relax and phonons

### relax
Create 1-optiDFTB folder and run 
```
relax(kPoints,geometry,temp,fmax)
```
* kPoints (list): number of kPoints for relaxation. default is [4,4,4]
* geometry: the structure file in the folder, ex.: "A.cif"
* temp (int): temperature (Defaults to 5 K)
* fmax (float): Maximum allowed force for convergence between atoms (Defaults to 0.01)

It will return optimized structure file "geo_end.gen" in the folder.
### phonons
