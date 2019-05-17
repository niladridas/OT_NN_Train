#!/bin/bash
#PBS -l nodes=1:ppn=10,mem=100gb,walltime=300:00:00
##
## Submit job to the no-limits queue (only available to members of "compnet" group)
#PBS -q high
##
## merge error and output to single file
#PBS -j oe
##
## Rename the output file
#PBS -o ex.o
##
## Rename the error file
##PBS -e example.sh.err
##
## Only run the job after the specified time, date is omittable as well as seconds
## however the hour and minut specifier need to be stated
##PBS -a [[[[CC]YY]MM]DD]hhmm.[SS]
##
## Specifies the arguments to pass to the script, valid options are like
##   qsub -F "myarg1 myarg2 myarg3=myarg3value"
##PBS -F "myarg1"
##
## send mail if the process aborts, when it begins, and
## when it ends (abe)
##PBS -m abe
##PBS -M <your username, e.g. sean.lawlor>
##
## Specifies custom environment variables to set for the script
##PBS -v TESTVAR="test value"
##
## Specifies the desired shell for this job
##PBS -S /bin/bash
##
## Job dependency: Specifies that a specific job must complete prior to this job.
## The following only shows if the job completed OK, otherwise see: FMI see:
## http://docs.adaptivecomputing.com/torque/5-0-1/Content/topics/torque/commands/qsub.htm
##PBS -W depend=afterok:<jobid>[:<jobid2>:<jobid3>:...]

######## Variables available to Torque jobs ##################
## PBS_JOBNAME	User specified jobname
## PBS_ARRAYID	Zero-based value of job array index for this job (in version 2.2.0 and later)
## PBS_GPUFILE	Line-delimited list of GPUs allocated to the job located in $TORQUE_HOME/aux/jobidgpu. Each line follows the following format:
##		<host>-gpu<number> For example, myhost-gpu1.
## PBS_O_WORKDIR	Job's submission directory
## PBS_TASKNUM	Number of tasks requested
## PBS_O_HOME	Home directory of submitting user
## PBS_MOMPORT	Active port for MOM daemon
## PBS_O_LOGNAME	Name of submitting user
## REM PBS_O_LANG	Language variable for job
## PBS_JOBCOOKIE	Job cookie
## PBS_JOBID	Unique pbs job id
## PBS_NODENUM	Node offset number
## PBS_NUM_NODES	Number of nodes allocated to the job
## PBS_NUM_PPN	Number of procs per node allocated to the job
## PBS_O_SHELL	Script shell
## PBS_O_HOST	Host on which job script is currently running
## PBS_QUEUE	Job queue
## PBS_NODEFILE	File containing line delimited list of nodes allocated to the job
## PBS_NP		Number of execution slots (cores) for the job
## PBS_O_PATH	Path variable used to locate executables within job script

# The CUDA device reserved for you by the batch system (just the integer ids)
# only needed if requesting gpus from above
CUDA_DEVICES=`cat $PBS_GPUFILE | rev | cut -d"-" -f1 | rev | cut -c "4" | sed ':a;N;$!ba;s/\n/ /g'`
# This is a "Matlab" friendly array of the device id's (starting at index 0, so you need to +1 to them in Matlab)
MATLAB_CUDA_ID_ARRAY="[$CUDA_DEVICES]"

# change your directory to where your code was submitted from
cd '/raid60/soumyasundar.pal/PFPFGMM';

## This function is an example of how to use MATLAB in a parallel/distributed
## way. i.e. parfor or spmd in MATLAB
function multipleMatlabWorkers {

# Start matlab with singleCompThread option so each matlab instance
# only uses 1 core
/usr/local/bin/matlab -nosplash -nodisplay -singleCompThread <<EOF

% Open the parallel pool with the number of cores requested
pool = parpool('local', $PBS_NUM_PPN);

run_all_model

% Delete the parallel pool
delete(gcp('nocreate'));

EOF
}

multipleMatlabWorkers
