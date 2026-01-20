#!/usr/bin/bash
#SBATCH --job-name=inPtest
#SBATCH --output=test_job.%j.out
#SBATCH --error=test_job.%j.err
#SBATCH --time=10:00:00
#SBATCH -p serc
#SBATCH -c 10
#SBATCH --mem=8GB

# change the working directory to current directory
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

# Write out some information on the job
echo Running on host `hostname`
echo Time is `date`
echo $SLURM_JOB_NAME 
### Define number of processors
echo This job has allocated $SLURM_NPROCS cpus

# Get the current and parent folder names
folder_name=$(basename "$PWD")
parent_folder=$(basename "$(dirname "$PWD")")

# Make output and backup directories
PDIR="/oak/stanford/groups/tchelepi/kimtaeho/Fervo/hbi"
ODIR="${PDIR}/${parent_folder}/${folder_name}"

if [ ! -d "$PDIR/${parent_folder}" ]; then
    mkdir -p "$PDIR/${parent_folder}"
fi

if [ ! -d "$ODIR" ]; then
    mkdir -p "$ODIR"
fi

echo "Output directory: $ODIR"

# if [ ! -e $IDIR ]; then
#     mkdir $IDIR 
# fi

rm -r output
ln -s $ODIR output
rm -r $ODIR/*

# rm -r input
# ln -s $IDIR input

# Tell me which nodes it is run on
echo " "
echo This jobs runs on the following processors:
echo $SLURM_JOB_NODELIST
echo " "

# Print out output directory
echo Output directory is $ODIR

ulimit -s unlimited
../lhbiem 3IstimTest.in
