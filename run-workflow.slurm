#!/bin/bash
## Slurm's SBATCH arguments (prefix with '#SBATCH '):
#SBATCH --job-name=beast_pype        # Job name
#SBATCH --output=beast_pype_%j.out    # Standard output file
#SBATCH --cpus-per-task=4             # Number of CPU cores per task
#SBATCH --time=0-01:00:00 # Maximum runtime (D-HH:MM:SS)
#SBATCH --mem=32G      # RAM to use.
### Depending on how slurm is configured on the HPC you are using you may have
### to add the name of a partition to use. If so remove the first # below and add
### the name after the -p:
##SBATCH -p    # Partition or queue name
##### More sbatch options can be found at https://slurm.schedmd.com/sbatch.html
###### The following are deactivated by the double hash instead of single.
#### I could not get the twi command below to work:
##SBATCH --mail-type=BEGIN,END,FAIL # Send email at job beginning, completion and failure.
##SBATCH --mail-user=  # Email address for notifications

#Load necessary conda environment using `source` NOT `conda`:
source activate beast_pype
workflow_notebook_path=$1
parameters_yaml=$2
job_out_file=beast_pype_${SLURM_JOB_ID}.out
parameters_string_yaml="$(cat $parameters_yaml)"
yaml() {
    python3 -c "import yaml;print(yaml.safe_load(open('$1'))$2)"
} # uses Pyyaml to extract variables from yaml.
overall_save_dir=$(yaml $parameters_yaml "['overall_save_dir']")
mkdir -p $overall_save_dir

if grep -Fxq 'specific_run_save_dir:' $parameters_yaml; then # Checks for phrase in file.
specific_run_save_dir=$(yaml $parameters_yaml "['specific_run_save_dir']")

else

specific_run_save_dir=$(date '+%Y-%m-%d_%H-%M-%S')
parameters_string_yaml=$(cat <<EOF
$parameters_string_yaml
specific_run_save_dir: $specific_run_save_dir
EOF
)

fi

sub_dir=$overall_save_dir/$specific_run_save_dir
mkdir -p $sub_dir
cache_dir="${sub_dir}/cache"
mkdir $cache_dir
echo "${SLURM_JOBID}" >> $cache_dir/slurm_job_ids.txt
workflow_file=`basename  $workflow_notebook_path`
papermill "$workflow_notebook_path" "$sub_dir/$workflow_file" -y "$parameters_string_yaml"

# Clean up
mv "$job_out_file" $sub_dir