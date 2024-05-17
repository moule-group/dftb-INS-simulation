# DFTB-INS-simulation

Use Density functional tight binding theory to simulate inelastic neutron scattering spectra.

## Installing phonopy
Create the environment in HPC (ex. NERSC) and install phonopy (I recommend numpy version=1.22.4 since there is bug if you have newer version of numpy using DFTB)  
```
conda create -p /global/common/software/xxxx/phonopy
conda install -c conda-forge phonopy numpy=1.22.4
```

## Installing DFTB 
Create the environment in your PC.
```
conda create -n dftb+
```
Or create the environment in HPC (ex. NERSC).
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
runDFTB(kPoints,sc,geometry=None,temp=5,fmax=1e-6,mode=None)
```
* kPoints (list): number of kPoints for relaxation. (Defaults to [8,8,8])
* sc (list): supercell size. (Defaults to [4,4,4])
* geometry: the structure file in the folder, ex.: "A.cif"
* temp (int): temperature (Defaults to 5 K)
* fmax (float): Maximum allowed force for convergence between atoms. (Defaults to 1e-6)

It will automatically find the structure file in the folder(*.cif, *.gen)

There are three modes, 
* relax 
* phonons
* force

### relax
```
python dftb.py -m relax
```
Create 1-optiDFTB folder and return optimized structure file "geo_end.gen" in the folder.
### phonons
```
python dftb.py -m phonons
```
The undistorted supercell is stored in geo.genS, while the required displacements are stored in files matching the pattern geo.genS-*

Create sub-folders to put geo.genS into the sub-folder.
```
parallel "mkdir {} && mv geo.genS-{} {}/geo_end.gen" ::: {001..xyz}
```
Copy dftb.py in 2-phonons folder & Run the single point energy (static) calculation by parallel command
```
parallel "cd {} && python ../dftb.py -m force -k k1 k2 k3" ::: {001..xyz}
```
### force

Reference: 
* https://phonopy.github.io/phonopy/dftb%2B.html
* https://dftbplus-recipes.readthedocs.io/en/latest/
