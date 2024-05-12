# DFTB-INS-simulation

Use Density functional tight binding theory to simulate inelastic neutron scattering spectra

## Installing DFTB 
Create the environment in your PC
```
conda create -n dftb
```
Or create the environment in supercomputer (ex. NERSC)
```
conda create -p /global/common/software/xxxx/dftb
```
Install DFTB by using mamba
```
conda activate dftb
mamba install 'dftbplus=*=nompi_*'
```
Download Slater-Koster files from DFTB website
```
* https://dftb.org/parameters/download/all-sk-files
```
