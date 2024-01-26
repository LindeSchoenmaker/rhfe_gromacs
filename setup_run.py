import pmx

from AZtutorial import AZtutorial

if __name__ == "__main__":
    print(pmx.__version__)
    # initialize the free energy environment object: it will store the main parameters for the calculations
    fe = AZtutorial( )

    # set the workpath
    fe.workPath = 'workpath_nosc'
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


    prepare = False
    if prepare:
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

    # set several parameters
    fe.JOBqueue = 'SLURM'
    fe.JOBsource = ['/etc/profile.d/modules.sh'] #,'/zfsdata/software/gromacs/2020.4/bin/GMXRC']
    fe.JOBmodules = ['gromacs/2022.1/gcc.8.4.0-cuda.11.7.1'] #['shared',' gmx_mpi','cuda11']
    fe.JOBexport = ['OMP_NUM_THREADS=8']
    fe.JOBgpu = True
    fe.JOBgmx = 'gmx mdrun'

    fe.prepare_jobscripts(simType='em')