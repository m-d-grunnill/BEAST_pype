#!/bin/bash
## Slurm's SBATCH arguments (prefix with '#SBATCH '):
#SBATCH --job-name=beast_pype        # Job name
#SBATCH --output=%x_%j.out  # Standard output file (the %x refers to --job-name and %j is jobid number)
#SBATCH --cpus-per-task=20             # Number of CPU cores per task
#SBATCH --time=0-01:00:00 # Maximum runtime (D-HH:MM:SS)
#SBATCH --mem=32G      # RAM to use.
### Depending on how slurm is configured on the HPC you are using you may have
### to add the name of a partition to use. If so remove the first # below and add
### the name after the -p:
#SBATCH --partition= # Partition name (if required by your HPC)
#SBATCH --account= # Account name for billing (if required by your HPC)
##### More sbatch options can be found at https://slurm.schedmd.com/sbatch.html
###### The following are deactivated by the double hash instead of single.
#### I could not get the this command below to work:
##SBATCH --mail-type=BEGIN,END,FAIL # Send email at job beginning, completion and failure.
##SBATCH --mail-user=  # Email address for notifications

source activate beast_pype #Load necessary conda environment using `source` NOT `conda`.
workflow=$1
parameters_yaml=$2
job_out_file=${SLURM_JOB_NAME}_${SLURM_JOB_ID}.out
yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
} # uses Pyyaml to extract variable_names from yaml.
overall_save_dir=$(yaml $parameters_yaml "['overall_save_dir']")
beast_pype run-workflow "$workflow" "$parameters_yaml"

# Clean up: move job_out_file to overall_save_dir if overall_save_dir exists.
if [ -d "$overall_save_dir" ]; then
  mv "$job_out_file" $overall_save_dir
fi