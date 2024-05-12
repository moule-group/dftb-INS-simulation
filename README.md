# DFTB-INS-simulation

Use Density functional tight binding theory to simulate inelastic neutron scattering spectra

## Installing DFTB 
Create the environment in your PC
```
conda create -n dftb+
```
Or create the environment in supercomputer (ex. NERSC)
```
conda create -p /global/common/software/xxxx/dftb+
```
Install DFTB by using mamba
```
conda activate dftb+
mamba install 'dftbplus=*=nompi_*'
```
Download Slater-Koster files from DFTB website
```
https://dftb.org/parameters/download/all-sk-files
```
If you want to access DFTB by ASE, you have to set environment variable in your .bashrc file
```
export DFTB_PREFIX="/path_to/slako/mio/mio-1-1/"
export ASE_DFTB_COMMAND="/path_to/dftb+/bin/dftb+ > PREFIX.out"
```
Add another environment variable in your .bashrc file in order to run dftb+ without activating envs everytime
```
export PATH"=/path_to/dftb+/bin:$PATH"
```
