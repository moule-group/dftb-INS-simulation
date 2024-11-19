###################### DFTB tool by ASE ###########################
# Create and modified by Ray started in May 2024 

# Updation 7/20/24: Adding DFT-D3 dispersion correction 
# Updation 7/23/24: 1. Adding Options_WriteChargesAsText='Yes', this can be used for partial charges in MD sim.; 2. Adding parser to set temperature value
# Updation 10/21/24 Adding DFTB-MD function
# Updation 11/5/24 Adding DFTB-MD dispersion function

###################################################################
import os
import argparse
import glob
from ase.io import read, write
import subprocess
import shutil
from ase.calculators.dftb import Dftb

# Constants

k_to_au = 0.316681534524639E-05 # 1 Kelvin
fs_to_au = 0.413413733365614E+02 # 1 fs
wn_to_au = 0.725163330219952E-06 # 1 cm-1

def phononsh(sc):
    """ Create phonon.sh file
    Args:
    sc (list): supercell size (Defaults to [2,2,2])
    #################################
    Return: phonon.sh file in the folder
    """
    s1 = int(sc[0])
    s2 = int(sc[1])
    s3 = int(sc[2])
    lines = ['#!/bin/bash\n'
             '\n',
             'mkdir 2-phonons\n',
             'cd 2-phonons\n',
             'cp ../1-relax/geo_end.gen .\n',
             'mv geo_end.gen geo.gen\n',
             f'phonopy -d --dim="{s1} {s2} {s3}" --dftb+\n',
             '\n',
             'exit 0']
    with open("phonon.sh", "w") as f:
        f.writelines(lines)
        
def forcesh(disp=False):
    """ Create force.sh file
    Args:
    disp (boolean): Dispersion correction using DFT-D3 (Defaults to False)
    #################################
    Return: force.sh file in the folder
    """
    if not disp:
        lines = ['#!/bin/bash\n'
             '\n',
             'if [ $# -eq 0 ]; then\n'
             '  echo "Usage: $0 <filename> $1 <folder number>"\n'
             '  exit 1\n'
             'fi'
             '\n'
             'cd $1\n'
             f'python ../../dftb.py -m force -k 1 1 1\n',
             '\n',
             'exit 0']
    else:
        lines = ['#!/bin/bash\n'
             '\n',
             'if [ $# -eq 0 ]; then\n'
             '  echo "Usage: $0 <filename> $1 <folder number>"\n'
             '  exit 1\n'
             'fi'
             '\n'
             'cd $1\n'
             f'python ../../dftb.py -m force -d -k 1 1 1\n',
             '\n',
             'exit 0']
    with open("force.sh", "w") as f:
        f.writelines(lines)
        
def restart_md():
    """ Create remd.sh file for restart md simulation
    Notes: If you want to continue md simulation, you can change ensemble from nvt to nve by modifying dftb_in.hsd in restart folder.
    This function now can only apply on same ensemble simulation (NVE)!
    Return: remd.sh file in the folder
    """
    lines = ['#!/bin/bash\n'
             '\n',
             'if [ $# -eq 0 ]; then\n'
             '  echo "Usage: $0 <filename> $1 <name of folder that simulation finish> $2 <name of the folder for next simulation>"\n'
             '  exit 1\n'
             'fi'
             '\n'
             'module load conda\n',
             'conda activate /global/common/software/m2734/ray/dftb+\n',
             'cd $1\n',
             'cp ../restart_filemaker.py .\n',
             'python restart_filemaker.py\n'
             'cd ..\n',
             'cp -r $1/restart .\n',
             'mv restart $2\n',
             'conda deactivate\n',
             '\n',
             'exit 0']
    with open("remd.sh", "w") as f:
        f.writelines(lines)
        
def runfile(name):
    """ Create run file for using NERSC HPC, this is written for GPU node.
    Args:
    name: The name of the run file.
    ##########################
    Return: run.sh file in the folder
    """
    
    with open(f'run.{name}', 'w') as f:
        cmds = ["#!/bin/bash\n",
                f"#SBATCH -J {name}\n",
                "#SBATCH -C gpu\n",
                "#SBATCH -A m2734\n",
                "#SBATCH -q regular\n",
                "#SBATCH -t 04:00:00\n",
                "#SBATCH -N 1\n",
                "#SBATCH -G 4\n",
                "#SBATCH --exclusive\n",
                "#SBATCH -o %x-%j.out\n",
                "#SBATCH -e %x-%j.err\n",
                "\n",
                "ulimit -s unlimited # for DFTB\n"
                "\n",
                "export OMP_NUM_THREADS=16\n",
                "export OMP_PLACES=threads\n",
                "export OMP_PROC_BIND=spread\n",
                "\n",
                f"srun -n 1 -c 128 --cpu_bind=cores -G 1 --gpu-bind=none dftb+ > {name}.out"]
        f.writelines(cmds)
        
def getGeometry(path):
    """ Using glob function in python to find the structure
    file in the current path
    The type of the structure files: "POSCAR", "*.gen"
    Args:
    path: The current directory (use os.getcwd())
    #########################################
    Return:
    file[0]: The structure file in the path
    """
    
    file = glob.glob(path + "/POSCAR") + glob.glob(path + "/*.gen") + glob.glob(path + "/*.xyz")
    
    return file[0]

def md(ensemble,geometry,kPoints,nsw,temp,disp=False):
    """ Run MD using DFTB calculatior 
    Args:
    ensemble (int): 1:NVT; 2:NVE
    geometry: The structure file in the folder, ex.: "POSCAR" or "*.gen"
    kPoints (list): Number of kPoints. (Defaults to [1,1,1])
    nsw (int): Number of MD steps
    temp (float): Temperature (Defaults to 150 K)
    disp (boolean): Dispersion correction using DFT-D3 (Defaults to False)
    """
    atoms = read(geometry)
    
    if not disp:
        calc_nvt = Dftb(label='nvt', 
                        Driver_='VelocityVerlet', # Verlet algorithm 
                        Driver_MovedAtoms='1:-1',
                        Driver_Steps=nsw,
                        Driver_MDRestartFrequency=4, # Info printed out frequency (every 4 steps)
                        Driver_Thermostat_='NoseHoover',
                        Driver_Thermostat_Temperature=temp*k_to_au,  # 1000K = 0.0031668 au. Here is 150K
                        Driver_Thermostat_CouplingStrength=3000*wn_to_au, # 3000 cm-1
                        Driver_TimeStep=0.25*fs_to_au, # 1fs = 41.3 a.u. Here is 0.25fs.
                        kpts=kPoints,
                        Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                        Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_MaxAngularMomentum_F='p',
                        Options_='',
                        Options_WriteChargesAsText='Yes',
                        Options_ReadChargesAsText='Yes',
                        )

        calc_nve = Dftb(label='nve', 
                        Driver_='VelocityVerlet', # Verlet algorithm 
                        Driver_MovedAtoms='1:-1',
                        Driver_Steps=nsw,
                        Driver_MDRestartFrequency=4, # Info printed out frequency (every 4 steps)
                        Driver_Thermostat_='None',
                        Driver_Thermostat_empty='',
                        Driver_TimeStep=0.25*fs_to_au, # 1fs = 41.3 a.u. Here is 0.25fs.
                        kpts=kPoints,
                        Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                        Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_MaxAngularMomentum_F='p',
                        Options_='',
                        Options_WriteChargesAsText='Yes',
                        Options_ReadChargesAsText='Yes',
                        )

    else: 
        calc_nvt = Dftb(label='nvt', 
                        Driver_='VelocityVerlet', # Verlet algorithm 
                        Driver_MovedAtoms='1:-1',
                        Driver_Steps=nsw,
                        Driver_MDRestartFrequency=4, # Info printed out frequency (every 4 steps)
                        Driver_Thermostat_='NoseHoover',
                        Driver_Thermostat_Temperature=temp*k_to_au,  # 1000K = 0.0031668 au. Here is 150K
                        Driver_Thermostat_CouplingStrength=3000*wn_to_au, # 3000 cm-1
                        Driver_TimeStep=0.25*fs_to_au, # 1fs = 41.3 a.u. Here is 0.25fs.
                        kpts=kPoints,
                        Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                        Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_MaxAngularMomentum_F='p',
                        Hamiltonian_Dispersion_='DftD3',
                        Hamiltonian_Dispersion_Damping='BeckeJohnson {',
                        Hamiltonian_Dispersion_Damping_a1='0.5719',
                        Hamiltonian_Dispersion_Damping_a2='3.6017',
                        Hamiltonian_Dispersion_s6='1.0',
                        Hamiltonian_Dispersion_s8='0.5883',
                        Options_='',
                        Options_WriteChargesAsText='Yes',
                        Options_ReadChargesAsText='Yes',
                        )
        
        calc_nve = Dftb(label='nve', 
                        Driver_='VelocityVerlet', # Verlet algorithm 
                        Driver_MovedAtoms='1:-1',
                        Driver_Steps=nsw,
                        Driver_MDRestartFrequency=4, # Info printed out frequency (every 4 steps)
                        Driver_Thermostat_='None',
                        Driver_Thermostat_empty='',
                        Driver_TimeStep=0.25*fs_to_au, # 1fs = 41.3 a.u. Here is 0.25fs.
                        kpts=kPoints,
                        Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                        Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                        Hamiltonian_MaxAngularMomentum_='',
                        Hamiltonian_MaxAngularMomentum_C='p',
                        Hamiltonian_MaxAngularMomentum_O='p',
                        Hamiltonian_MaxAngularMomentum_H='s',
                        Hamiltonian_MaxAngularMomentum_N='p',
                        Hamiltonian_MaxAngularMomentum_S='d',
                        Hamiltonian_MaxAngularMomentum_F='p',
                        Hamiltonian_Dispersion_='DftD3',
                        Hamiltonian_Dispersion_Damping='BeckeJohnson {',
                        Hamiltonian_Dispersion_Damping_a1='0.5719',
                        Hamiltonian_Dispersion_Damping_a2='3.6017',
                        Hamiltonian_Dispersion_s6='1.0',
                        Hamiltonian_Dispersion_s8='0.5883',
                        Options_='',
                        Options_WriteChargesAsText='Yes',
                        Options_ReadChargesAsText='Yes',
                        )
        
    if ensemble == 1:
        atoms.calc = calc_nvt
        calc_nvt.write_input(atoms) 
     
    elif ensemble == 2:
        atoms.calc = calc_nve
        calc_nve.write_input(atoms) 

def relax(kPoints,geometry,temp,fmax,disp=False):
    """ Relax structure by using DFTB
    Documentation: 
    Driver: determines how the geometry should be changed 
    during the calculation.
    Driver_Moveatoms="1:-1": Move all atoms in the system
    Hamiltonian_MaxAngularMomentum_C='p':
    the value of the highest orbital angular momentum each element
    ##############################################
    Args:
    kPoints (list): number of kPoints for relaxation. (Defaults to [1,1,1], it should be increased.)
    geometry: the structure file in the folder, ex.: "POSCAR" or "*.gen"
    temp (float): temperature (Defaults to 150 K)
    fmax (float): Maximum allowed force for convergence between atoms (Defaults to 5e-4)
    disp (boolean): Dispersion correction using DFT-D3 (Defaults to False)
    """
    atoms = read(geometry)
    
    if not disp:
        calc = Dftb(label='relax', 
                Driver_='ConjugateGradient', # relax algorithm (This command will be change in later version)
                Driver_MovedAtoms='1:-1',
                Driver_MaxForceComponent=fmax, # converge if forces on all atoms < fmax
                Driver_MaxSteps=1000,
                Driver_LatticeOpt='Yes',
                Driver_Isotropic='Yes',
                kpts=kPoints,
                Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force < tolerance
                Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}'.format(T=temp), # Set this temp = 300K
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_F='p',
                Options_='',
                Options_WriteChargesAsText='Yes'
                )
    
    else:
        calc = Dftb(label='relax', 
                Driver_='ConjugateGradient', # relax algorithm (This command will be change in later version)
                Driver_MovedAtoms='1:-1',
                Driver_MaxForceComponent=fmax, # converge if forces on all atoms < fmax
                Driver_MaxSteps=1000,
                Driver_LatticeOpt='Yes',
                Driver_Isotropic='Yes',
                kpts=kPoints,
                Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force below tolerance
                Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}'.format(T=temp), # Set this temp = 300K
                Hamiltonian_MaxAngularMomentum_='', 
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_F='p',
                Hamiltonian_Dispersion_='DftD3',
                Hamiltonian_Dispersion_Damping='BeckeJohnson {',
                Hamiltonian_Dispersion_Damping_a1='0.5719',
                Hamiltonian_Dispersion_Damping_a2='3.6017',
                Hamiltonian_Dispersion_s6='1.0',
                Hamiltonian_Dispersion_s8='0.5883',
                Options_='',
                Options_WriteChargesAsText='Yes'
                )

    #atoms.set_calculator(calc) 
    #calc.calculate(atoms)
    #with open('dftb_in.hsd', 'w') as f:
    #calc.write_dftb_in(f) 
    calc.write_input(atoms) 
    
def singleE(geometry,kPoints,temp,disp=False):
    """ Single point energy calculation and extract force by using DFTB,
    this step is for getting phonons with phonopy.
    Args:    
    geometry: the structure file in the folder, ex.: "POSCAR or *.gen"
    kPoints (list): kPoints list for relax. default is [1,1,1]
    temp (float): temperature (Defaults to 150 K)
    disp (boolean): Dispersion correction using DFT-D3 (Defaults to False)
    """
    atoms = read(geometry)
    
    if not disp:
        calc = Dftb(label='singleE', 
                kpts=kPoints,
                Analysis_='',
                Analysis_CalculateForces='Yes', # This command changed in dftb+ 24.1
                Options_="",
                Options_WriteResultsTag='Yes',
                Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force below tolerance
                Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}'.format(T=temp),
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_F='p')
    else:
        calc = Dftb(label='singleE', 
                kpts=kPoints,
                Analysis_='',
                Analysis_CalculateForces='Yes', # This command changed in dftb+ 24.1
                Options_='',
                Options_WriteResultsTag='Yes',
                Hamiltonian_SCC='Yes', # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-6, # Stop if maximal force below tolerance
                Hamiltonian_Filling='Fermi {{Temperature [Kelvin] = {T} }}'.format(T=temp),
                Hamiltonian_MaxAngularMomentum_='',
                Hamiltonian_MaxAngularMomentum_C='p',
                Hamiltonian_MaxAngularMomentum_O='p',
                Hamiltonian_MaxAngularMomentum_H='s',
                Hamiltonian_MaxAngularMomentum_N='p',
                Hamiltonian_MaxAngularMomentum_S='d',
                Hamiltonian_MaxAngularMomentum_F='p',
                Hamiltonian_Dispersion_='DftD3',
                Hamiltonian_Dispersion_Damping='BeckeJohnson {',
                Hamiltonian_Dispersion_Damping_a1='0.5719',
                Hamiltonian_Dispersion_Damping_a2='3.6017',
                Hamiltonian_Dispersion_s6='1.0',
                Hamiltonian_Dispersion_s8='0.5883')
        
    with open('dftb_in.hsd', 'w') as f:
        calc.write_dftb_in(f) 
    os.system("dftb+ 1>> forces.out 2>> forces.err") # Run this in order to get result.tag
    
def runDFTB(ensemble,fmax,kPoints,mode,nsw,sc,temp,disp=False):
    """ Run DFTB, there are 5 modes 1. relax 2. phonons. Need to specify
    the mode to run this function. For 3. force, it will run by command "bash force.sh". 4. md. 5. restart md.
    Args:
    ensemble: if mode = 4, need specify the ensemble (NVT or NVE).
    fmax: Defaults to 5e-4
    kPoints: Defaults to [1,1,1]
    mode (int): Need to specify
    nsw (int): Number of steps for MD simulation
    sc (list): Defaults to [2,2,2]
    temp (float): Defaults to 150K
    disp (boolean): Dispersion correction using DFT-D3 (Defaults to False)
    """
    path = os.getcwd()
    
    if mode == 1:
        geometry = getGeometry(path)
        os.mkdir(path + "/1-relax")
        os.chdir(os.path.join(path, "1-relax")) # Change directory to 1-relax
        try: 
            relax(kPoints,geometry,temp,fmax,disp)
            runfile(name='relax')
        finally:
            os.chdir(path)
            
    if mode == 2:
        phononsh(sc=sc)
        subprocess.run(["bash" ,"phonon.sh"])
        os.chdir(os.path.join(path, "2-phonons"))
        forcesh(disp)
    
    if mode == 3: 
        geometry = getGeometry(path)
        singleE(kPoints,geometry,temp,disp)
        
    if mode == 4:
        geometry = getGeometry(path)
        
        if ensemble == 1:
            os.mkdir(path + "/1-nvtmd")
            os.chdir(os.path.join(path, "1-nvtmd")) 
            try:     
                md(ensemble,geometry,kPoints,nsw,temp,disp)
                runfile(name="nvtmd")
            finally:
                os.chdir(path)
                
        elif ensemble == 2:
            os.mkdir(path + "/1-nvemd")
            os.chdir(os.path.join(path, "1-nvemd")) 
            try:     
                md(ensemble,geometry,kPoints,nsw,temp,disp)
                runfile(name="nvemd")
            finally:
                os.chdir(path)
            
    if mode == 5:
        restart_md()
        runfile(name="restart")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser() # Create the parser
    parser.add_argument("-m", "--mode", type=int, help="1:relax, 2:phonons, 3:force, 4:md, 5:re_md") # Add an argument: mode
    parser.add_argument("-f", "--fmax", type=float, default=5e-4, help="converge if forces on all atoms < fmax") # Add an argument: force
    parser.add_argument("-t", "--temp", type=float, default=150.0, help="temperature in Kelvin") # Add an argument: temp
    parser.add_argument("-k", "--kpoints", type=int, nargs=3, default=[1,1,1], help="k1, k2, k3") # Add an argument: kPoints
    parser.add_argument("-s", "--supercell", type=int, nargs=3, default=[2,2,2], help="supercell size") # Add an argument: supercell   
    parser.add_argument("-d", "--dispersion", action='store_true', help="Dispersion correction (DFT-D3)") # Add an argument: dispersion 
    parser.add_argument("-e", "--ensemble", type=int, default=1, help="1:NVT or 2:NVE") # Add an argument: ensemble
    parser.add_argument("-n", "--nsw", type=int, default=4000, help="Total MD steps") # Add an argument: nsw
    args = parser.parse_args() # Parse the argument
    runDFTB(ensemble=args.ensemble,fmax=args.fmax,kPoints=args.kpoints,mode=args.mode,nsw=args.nsw,sc=args.supercell,temp=args.temp,disp=args.dispersion)
 
