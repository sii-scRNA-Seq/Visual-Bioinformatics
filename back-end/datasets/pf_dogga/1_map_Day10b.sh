#!/bin/bash -l

############# SLURM SETTINGS #############
#SBATCH --account=project0001   # account name (mandatory), if the job runs under a project then it'll be the project name, if not then it should =none
#SBATCH --job-name=map_Day10b        # some descriptive job name of your choice
#SBATCH --output=%x-%j.out      # output file name will contain job name + job ID
#SBATCH --error=%x-%j.err       # error file name will contain job name + job ID
#SBATCH --partition=nodes        # which partition to use, default on MARS is “nodes"
#SBATCH --time=0-12:00:00       # time limit for the whole run, in the form of d-hh:mm:ss, also accepts mm, mm:ss, hh:mm:ss, d-hh, d-hh:mm
#SBATCH --mem=64G                # memory required per node, in the form of [num][M|G|T]
#SBATCH --nodes=1               # number of nodes to allocate, default is 1
#SBATCH --ntasks=1              # number of Slurm tasks to be launched, increase for multi-process runs ex. MPI
#SBATCH --cpus-per-task=32       # number of processor cores to be assigned for each task, default is 1, increase for multi-threaded runs
#SBATCH --ntasks-per-node=1     # number of tasks to be launched on each allocated node

cd /mnt/data/project0001/2117532m/SCAMPI_PF_dataset
module load apps/cellranger/7.2.0

cellranger count --disable-ui --fastqs=/mnt/data/project0001/2117532m/SCAMPI_PF_dataset/MCA_data/ERR11471992,/mnt/data/project0001/2117532m/SCAMPI_PF_dataset/MCA_data/ERR11471993,/mnt/data/project0001/2117532m/SCAMPI_PF_dataset/MCA_data/ERR11471994,/mnt/data/project0001/2117532m/SCAMPI_PF_dataset/MCA_data/ERR11471995 --transcriptome=PF --sample=ERR11471992,ERR11471993,ERR11471994,ERR11471995 --id=MCA_PF_counts_Day10b --localmem 64 --localcores 32 
