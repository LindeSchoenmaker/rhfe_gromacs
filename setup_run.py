import argparse

import pmx

from AZtutorial import AZtutorial

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", help="which files to produce",
                type=str, choices=['initial', 'em', 'equil_nvt', 'equil_npt', 'production'])
parser.add_argument("-wp", "--workPath", help="the workpath",
                type=str, default="workpath")
parser.add_argument("-p", "--JOBpartition", help="which partition to use", default ='free-gpu',
                type=str, choices=['free', 'free-gpu', 'standard', 'gpu'])
parser.add_argument("-t", "--JOBsimtime", help="simulation time in hours", default =3,
                type=int)

if __name__ == "__main__":
    print(pmx.__version__)

    args = parser.parse_args()
    # initialize the free energy environment object: it will store the main parameters for the calculations
    fe = AZtutorial( )

    # set the workpath
    fe.workPath = args.workPath
    # set the path to the molecular dynamics parameter files
    fe.mdpPath = 'input/mdppath/files'
    # set the number of replicas (several repetitions of calculation are useful to obtain reliable statistics)
    fe.replicas = 1
    # provide the path to the protein structure and topology
    fe.proteinPath = 'None'
    # provide the path to the folder with ligand structures and topologies
    fe.ligandPath = 'input/ligands'
    # provide edges
    fe.edges = [ ['to_', 'int'] ] #, ['int','ref'], ['to_', 'ref'] ]
    # finally, let's prepare the overall free energy calculation directory structure
    fe.prepareFreeEnergyDir( )

    # set several parameters
    fe.JOBqueue = 'SLURM'
    fe.JOBsource = ['/etc/profile.d/modules.sh'] #,'/zfsdata/software/gromacs/2020.4/bin/GMXRC']
    fe.JOBmodules = ['gromacs/2022.1/gcc.8.4.0-cuda.11.7.1'] #['shared',' gmx_mpi','cuda11']
    fe.JOBexport = ['OMP_NUM_THREADS=8']
    fe.JOBgpu = True
    fe.JOBgmx = 'gmx mdrun'
    fe.JOBpartition = args.JOBpartition

    fe.JOBsimtime = args.JOBsimtime

    if args.output == 'initial':
        # this command will map the atoms of all edges found in the 'fe' object
        # bVerbose flag prints the output of the command
        fe.atom_mapping(bVerbose=False)
        #construct hybrid topology
        fe.hybrid_structure_topology(bVerbose=False, bSeparate=True)
        #assemble ligand+water systems
        fe.assemble_systems( )
        #build box, solvate
        fe.boxWaterIons( )

        #prepare simulation
        fe.prepare_simulation( simType='em' )
        fe.prepare_jobscripts(simType='em')
    else:
        fe.prepare_simulation( simType=args.output )
        fe.prepare_jobscripts(simType=args.output)