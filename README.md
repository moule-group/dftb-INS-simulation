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
runDFTB(ensemble,fmax,kPoints,mode,nsw,sc,temp,disp=False)
```
* ensemble: if mode = 4, need specify the ensemble (NVT or NVE).
* fmax (float): Maximum allowed force for convergence between atoms. (Defaults to 5e-4)
* kPoints (list): number of kPoints for relaxation. (Defaults to [1,1,1], it should be increased for accuracy)
* mode (int): Need to specify
* nsw (int): Number of steps for MD simulation
* sc (list): supercell size. (Defaults to [2,2,2])
* temp (float): Defaults to 150K
* disp (boolean): Dispersion correction using DFT-D3 (Defaults to False if not specifiy in command)

It will automatically find the structure file in the folder(POSCAR, /*.gen)

There are 5 modes, 
* 1-relax: Relax crystal structures
* 2-phonons: Create finite displacement structure files and put into separate folders.
* 3-force: Single point energy calculation and extract force by using DFTB, this step is for getting phonons with phonopy.
* 4-md: DFTB-MD simulation
* 5-restart_md: Restart DFTB-MD simulation if necessary
 
### relax
```
python dftb.py -m relax
```
Create 1-optiDFTB folder and return optimized structure file "geo_end.gen" in the folder.
### phonons
```
python dftb.py -m phonons
```
Write phonon.sh file and run this file. Return the undistorted supercell is stored in geo.genS, while the required displacements are stored in files matching the pattern geo.genS-*

Create sub-folders to put geo.genS into the sub-folder (using parallel).
```
parallel "mkdir {} && mv geo.genS-{} {}/geo_end.gen" ::: {001..xyz}
```

### force
Copy dftb.py in 2-phonons folder & Run the single point energy (static) calculation by parallel command
```
parallel --delay .2 bash force.sh ::: {001..xyz}
```
After finishing all the force calculation, run the command below
```
phonopy -f {001..xyz}/results.tag --dftb+
```
Create mesh.conf
```
ATOM_NAME = C H
DIM = 2 2 2
MP = 8 8 8
EIGENVECTORS = .TRUE.
   ```
Calculate phonons and DOS. We will get "phonopy.yaml" and "mesh.yaml" file.
```
phonopy -p mesh.conf
```

#### md
For NVT ensemble, here we use Nose-Hoover thermostat, set -e 1.
```
python dftb.py -m 4 -e 1
```
For NVE ensemble, set -e 2.
```
python dftb.py -m 4 -e 2
```

#### restart_md

#### appendix

srun -n --ntasks: number of tasks, equals to ntasks_per_node x num_of_node

Total Logical Processors = Number of Physical Cores × Threads per Core; (in Perlmutter, each node has 2 64 core AMD CPUs and each AMD CPU has 2 hardware threads)


srun -c --cpus-per-task: equals to num_of_cpu (2x2x64) / ntasks_per_node

Reference: 
* https://phonopy.github.io/phonopy/dftb%2B.html
* https://dftbplus-recipes.readthedocs.io/en/latest/
