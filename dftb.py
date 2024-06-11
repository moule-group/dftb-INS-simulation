# Create and modified by Ray started in May 2024 (last updated 5/28/24)
import os
import argparse
#import oclimax
import glob
from ase.io import read
import subprocess
#from ase.optimize import BFGS, LBFGS
from ase.calculators.dftb import Dftb

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
        
def forcesh():
    """ Create force.sh file
    #################################
    Return: force.sh file in the folder
    """
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
    with open("force.sh", "w") as f:
        f.writelines(lines)
        
def getGeometry(path):
    """ Using glob function in python to find the structure
    file in the current path
    The type of the structure files: "POSCAR", "gen"
    Args:
    path: The current directory (use os.getcwd())
    #########################################
    Return:
    file: The structure file in the path
    """
    
    file = glob.glob(path + "/POSCAR") + glob.glob(path + "/*.gen")
    
    return file[0]

def relax(kPoints,geometry,temp,fmax):
    """ Relax structure by using DFTB
    Documentation: 
    Driver: determines how the geometry should be changed 
    during the calculation.
    Driver_Moveatoms="1:-1": Move all atoms in the system
    Hamiltonian_MaxAngularMomentum_C='p':
    the value of the highest orbital angular momentum each element
    ##############################################
    Args:
    kPoints (list): number of kPoints for relaxation. (Defaults to [8,8,8])
    geometry: the structure file in the folder, ex.: "A.cif"
    temp (int): temperature (Defaults to 5 K)
    fmax (float): Maximum allowed force for convergence between atoms (Defaults to 1e-3)
    """
    
    atoms = read(geometry)
    calc = Dftb(label="relax", 
                Driver_="ConjugateGradient", # relax algorithm (This command will be change in later version)
                Driver_MovedAtoms="1:-1",
                Driver_MaxForceComponent=fmax, # converge if forces on all atoms < fmax
                Driver_MaxSteps=1000,
                Driver_LatticeOpt="Yes",
                Driver_Isotropic="Yes",
                kpts=kPoints,
                Hamiltonian_SCC="Yes", # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-8, # Stop if maximal force < tolerance
                Hamiltonian_Filling="Fermi {{Temperature [Kelvin] = {T} }}".format(T=temp),
                Hamiltonian_MaxAngularMomentum_="",
                Hamiltonian_MaxAngularMomentum_C="p",
                Hamiltonian_MaxAngularMomentum_O="p",
                Hamiltonian_MaxAngularMomentum_H="s",
                Hamiltonian_MaxAngularMomentum_N="p",
                Hamiltonian_MaxAngularMomentum_S="d",
                Hamiltonian_MaxAngularMomentum_F="p")

    atoms.set_calculator(calc) 
    calc.calculate(atoms)
    
def singleE(kPoints,geometry,temp):
    """ single point energy calculation and extract force by using DFTB,
    this step is for getting phonons with phonopy.
    Args:
    kPoints (list): kPoints list for relax. default is [8,8,8] (It can be reduced!)
    geometry: the structure file in the folder, ex.: "A.cif"
    temp (int): temperature default is 5 units: [K]
    """
    atoms = read(geometry)
    calc = Dftb(label="singleE", 
                kpts=kPoints,
                Analysis_="",
                Analysis_CalculateForces="Yes", # This command changed in dftb+ 24.1
                Options_="",
                Options_WriteResultsTag="Yes",
                Hamiltonian_SCC="Yes", # self-consistent redistribution of Mulliken charges
                Hamiltonian_SCCTolerance=1e-8, # Stop if maximal force below tolerance
                Hamiltonian_Filling="Fermi {{Temperature [Kelvin] = {T} }}".format(T=temp),
                Hamiltonian_MaxAngularMomentum_="",
                Hamiltonian_MaxAngularMomentum_C="p",
                Hamiltonian_MaxAngularMomentum_O="p",
                Hamiltonian_MaxAngularMomentum_H="s",
                Hamiltonian_MaxAngularMomentum_N="p",
                Hamiltonian_MaxAngularMomentum_S="d",
                Hamiltonian_MaxAngularMomentum_F="p")
    
    with open('dftb_in.hsd', 'w') as f:
        calc.write_dftb_in(f) 
    os.system("dftb+ 1>> forces.out 2>> forces.err") # Run this in order to get result.tag
    
    
def runDFTB(kPoints=[4,4,4],sc=[2,2,2],geometry=None,temp=5,fmax=1e-3,mode=None):
    """ run DFTB, there are two modes 1. relax 2. phonons. Need to specify
    the mode to run this function
    """
    path = os.getcwd()
    
    if mode == "relax":
        geometry = getGeometry(path)
        os.mkdir(path + "/1-relax")
        os.chdir(os.path.join(path, "1-relax")) # Change directory to 1-relax
        try: 
            relax(kPoints=kPoints,geometry=geometry,temp=temp,fmax=fmax)
        finally:
            os.chdir(path)
            
    elif mode == "phonons":
        phononsh(sc=sc)
        subprocess.run(["bash" ,"phonon.sh"])
        os.chdir(os.path.join(path, "2-phonons"))
        forcesh()
    
    elif mode == "force": 
        geometry = getGeometry(path)
        singleE(kPoints=kPoints,geometry=geometry,temp=temp)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser() # Create the parser
    parser.add_argument("-m", "--mode", help="relax, phonons, force") # Add an argument: mode
    parser.add_argument("-f", "--fmax", type=float, default=1e-3, help="converge if forces on all atoms < fmax") # Add an argument: force
    parser.add_argument("-k", "--kpoints", type=int, nargs=3, default=[4,4,4], help="k1, k2, k3") # Add an argument: kPoints
    parser.add_argument("-s", "--supercell", type=int, nargs=3, default=[2,2,2], help="supercell size") # Add an argument: supercell       
    args = parser.parse_args() # Parse the argument
    runDFTB(kPoints=args.kpoints, mode=args.mode, sc=args.supercell, fmax=args.fmax)
