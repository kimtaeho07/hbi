#!/usr/bin/bash
#SBATCH --job-name=test_job
#SBATCH --output=test_job.%j.out
#SBATCH --error=test_job.%j.err
#SBATCH --time=05:00
#SBATCH -p normal
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

### Make output and back up directories 
ODIR=/oak/stanford/groups/tchelepi/kimtaeho/Fervo/hbi/hbi_gaussianSource/pSig1
if [ ! -e $ODIR ]; then
    mkdir $ODIR 
fi
ln -s $ODIR output
rm -r $ODIR/*

# Tell me which nodes it is run on
echo " "
echo This jobs runs on the following processors:
echo $SLURM_JOB_NODELIST
echo " "

# Print out output directory
echo Output directory is $ODIR

../lhbiem pSig1.in
