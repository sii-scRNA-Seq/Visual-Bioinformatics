#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0001   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=makereference        # some descriptive job name of your choice
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=0-24:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=8G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=48       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

############# LOADING MODULES (optional) #############
#module load apps/xxx
#module load libs/xxx

cd /mnt/data/project0001/2117532m/SCAMPI_PF_dataset

#source cellranger-7.2.0/sourceme.bash
module load apps/agat/1.2.0

############# MY CODE #############
echo "Hello from $SLURM_JOB_NODELIST"
echo "PWD: $PWD"

echo "Converting GFF to GTF"
agat_convert_sp_gff2gtf.pl --gff PlasmoDB-68_Pfalciparum3D7.gff --output PlasmoDB-68_Pfalciparum3D7.gtf

rm *.agat.log
module unload apps/agat/1.2.0
echo "Making Reference"
module load apps/cellranger/7.2.0
cellranger mkref --genome=PF --fasta=PlasmoDB-68_Pfalciparum3D7_Genome.fasta --genes=PlasmoDB-68_Pfalciparum3D7.gtf --nthreads=48
