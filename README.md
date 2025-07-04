# BEAST_pype: A pipeline for phylodynamics using BEAST2.

## IMPORTANT NOTES
* With this early access version of beast_pype the beast_pype python package is
  installed with the `-e` flag. This means it is intsalled in an editable development mode.
  **MEANING IF YOU MAKE ANY CHANGES TO THE FILES IN  `src/beast_pype` THIS WILL
  AFFECT YOUR beast_pype INSTALLATION**
* Currently, beast_pype is reliant on [slurm](https://slurm.schedmd.com/overview.html). 
If you are running a beast_pype workflow on a machine without [slurm](https://slurm.schedmd.com/overview.html) it can be installed  via (conda)[https://anaconda.org/conda-forge/slurm]. However, the 
initial configuration of slurm may be quite involved.
* For licensing reasons beast_pype uses the version of BEAST 2 that is available via [bioconda](https://anaconda.org/bioconda/beast2), 2.6 as 2025-Jun-25. Template beast2 xmls from other versions of may not work. BEAST 2.7.7 is available on the conda channel [millerjeremya](https://anaconda.org/millerjeremya/beast2). However,
as of 2025-06-026 I can't tell if this channel is under an open license (like conda-forge or bioconda). Some conda channels charge organisation of >200 staff for usage (https://www.datacamp.com/blog/navigating-anaconda-licensing). 

## Initial Setup:
You should only need to carry out the following steps the first time you run the pipeline 
1. Clone the beast_pype repo: 
  ```bash
    git clone https://gitlab.cscscience.ca/siage/projects/BEAST_pipe.git
 ```
2. Install the beast_pype conda environment via conda or mamba:
```bash
conda env create -f requirements.yml
```
```bash
mamba env create -f requirements.yml
```
3. For Jupyter you will need to add a beast_pype environment **python** kernel and a **bash** kernel:
```bash
conda activate beast_pype
python -m ipykernel install --user --name=beast_pype
python -m bash_kernel.install --user
```
4. The beast_pype package itself will need to be installed:
```bash
conda activate beast_pype
cd LOCATION_YOU_CLONED_THE_BEAST_PYPE_REPO_TO
pip install -e .
```
5. You will need to install BEAST 2 packages for the BEAST that a template BEAST 2 xml uses (see )
    The command below demonstrates how to install the BDSKY package used in the example:
    ```
    conda activate beast_pype
    packagemanager -add BDSKY
    ```
6. For using the widgets (interactive GUI) in Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb notebook you may need to run the commands below:
    ```
    conda activate beast_pype
    jupyter nbextension enable --py --sys-prefix widgetsnbextension
    jupyter labextension install @jupyter-widgets/jupyterlab-manager
    ```


## Launching beast_pype Workflows via Slurm's sbatch 

Beast_pype workflows are launched using Slurm's sbatch via the following form:
```bash
sbatch run-workflow.slurm PATH_TO_WORKFLOW.ipnyb PATH_TO_PARAMETERS.yml
```
`run-workflow.slurm` is a specialised [shell script for running sbatch jobs](https://hpc-uit.readthedocs.io/en/latest/jobs/examples.html). The
top of `run-workflow.slurm` has `#SBATCH` arguments that you may
wish to alter or add to (see  [slurm's sbatch documentation](https://slurm.schedmd.com/sbatch.html)) depending on a particular run of a workflow.

**BEFORE RUNNING THE EXAMPLES BELOW** from the administrator of the HPC you
are using you will need to learn the name of the partition to use for your slurm jobs. 
Once you have acquired this you will need to alter line 5 of `run-workflow.slurm`  to read:
```bash
#SBATCH -p NAME_OF_PARTITION_TO_USE
```
Likewise, any parameter ymls used for running a workflow (such as
`parameters/Test-BDSKY-serial.yml` and `parameters/Test-Generic.yml`) need to 
have the line `partition: NAME_OF_PARTITION_TO_USE` added.

### Example Run of BDSKY-serial Workflow

This workflow is intended for use with the Birth-Death Skyline (BDSKY)-serial model (see https://www.beast2.org/files/bdskytutorialv2.1.1.pdf). This will only work with a template
BDSky-serial xml or providing a ready-to-go BDSky-serial xml to use. A template
BDSky-serial xml has the sequences from a fasta file inserted (depending on other
settings an initial tree added, the Rt dimensions and origin prior altered) before being run by BEAST 2.  A ready-to-go BDSky-serial xml will simply be run by BEAST 2. See parameters section at the top of [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb) for more detailed
information on workflow parameters.  

The [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb) along with other workflows
can be found in the `workflows` directory. To run this from commandline via slurms `sbatch`
enter the following command
```bash
cd LOCATION_YOU_CLONED_THE_BEAST_PYPE_REPO_TO
sbatch run-workflow.slurm workflows/BDSKY-serial.ipynb parameters/Test-BDSKY-serial.yml
```
This will run the workflow with the parameters set in `parameters/Test-BDSKY-serial.yml`. 
You can alter them for usage or alter a copy of `parameters/Test-BDSKY-serial.yml`.
A description of the parameters used in [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb)
is at the top of the notebook. **NOTE** the default values for `chain_length`,
`store_state_every`, `trace_log_every`, `tree_log_every`, `screen_log_every` and `store_state_every` are extremely low. These low values are intended for testing 
that a beast_pype workflow works on your set-up (taking a few minutes to run). **IT 
IS RECOMMEDED THAT** you up the values of `chain_length` (>=10,000,000),
`store_state_every` (>=1000), `trace_log_every` (>=1000), `tree_log_every` (>=1000), `screen_log_every` (>=1000) and `store_state_every`(>=1000) when properly running this workflow. 

Running the above command will create the directory `example_runs_of_BSDKY`
containing a **time stamped** folder Shortly after the BDSKY-serial workflow is launched
a running copy of [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb) should appear.
The top code cell of this notebook will have a new code cell inserted below which replaces the default parameter values
in [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb) with the  values specified in
`parameters/Test-BDSKY-serial.yml`.  The 'time stamped' folder will gradually fill with running
copies of the notebooks (in the order):
1. Phase-2i-IQTree.ipynb
2. Phase-2ii-TreeTime-and-Down-Sampling.ipynb
3. Phase-3-Gen-BDSKY-xml.ipynb
4. Phase-4-Running-BEAST.ipynb

Lastly, `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb` should appear in the
**time stamped** folder  Once all the MCMC BEAST 2 chains have finished running the copy of `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb` needs to be done manually. 
The progress of the MCMC BEAST 2 chains can be checked via openning the files `run-1.out`,
`run-2.out`, `run-3.out` and `run-4.out` or the command:
```bash
   squeue --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R" --me
```
The slurm jobs running MCMC BEAST 2 chains will be named `beast_full_sample_run_1`,
`beast_full_sample_run_2`, `beast_full_sample_run_3` and `beast_full_sample_run_4`.
Competing the manual run of the copy of `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb`
produces the output report notebook `BEAST_pype-Report.ipynb` and html version `BEAST_pype-Report.html`.

### Example Run of Genric Workflow

This workflow is intended for use with any beast 2 xml. This will work with any 
BEAST 2 xml template or ready-to-go BEAST 2 xml. A template
BEAST 2 xml has the sequences from a fasta file inserted, and depending on other
settings an initial tree added before being run by BEAST 2.  A 
ready-to-go BDSky-serial xml will simply be run by BEAST 2. See parameters section
at the top of [Generic workflow notebook](workflows/Generic.ipynb) for more detailed information on
workflow parameters.  

The [Generic workflow notebook](workflows/Generic.ipynb) along with other workflows
can be found in the `workflows` directory. To run this from commandline via slurms `sbatch`
enter the following commands 
```bash
cd LOCATION_YOU_CLONED_THE_BEAST_PYPE_REPO_TO
sbatch run-workflow.slurm workflows/Generic.ipynb parameters/Test-Generic.yml
```
This will run the workflow with the parameters set in `parameters/Test-Generic.yml`. 
You can alter them for usage or alter a copy of `parameters/Test-Generic.yml`.
A description of the parameters used in [Generic workflow notebook](workflows/Generic.ipynb)
is at the top of the notebook. **NOTE** the default values for `chain_length`,
`store_state_every`, `trace_log_every`, `tree_log_every`, `screen_log_every` and `store_state_every` are extremely low. These low values are intended for testing 
that a beast_pype workflow works on your set-up (taking a few minutes to run). **IT 
IS RECOMMEDED THAT** you up the values of `chain_length` (>=10,000,000),
`store_state_every` (>=1000), `trace_log_every` (>=1000), `tree_log_every` (>=1000), `screen_log_every` (>=1000) and `store_state_every`(>=1000) when properly running this workflow. 

Running the above command will create the directory `example_runs_of_BSDKY`
containing a **time stamped** folder Shortly after the Generic workflow is launched
a running copy of [Generic workflow notebook](workflows/Generic.ipynb) should appear.
The top code cell of this notebook will have a new code cell inserted below which replaces the default parameter values
in [Generic workflow notebook](workflows/Generic.ipynb) with the  values specified in
`parameters/Test-Generic.yml`.  The 'time stamped' folder will gradually fill with running
copies of the notebooks (in the order):
1. Phase-2i-IQTree.ipynb
2. Phase-2ii-TreeTime-and-Down-Sampling.ipynb
3. Phase-3-Gen-Generic-xml.ipynb
4. Phase-4-Running-BEAST.ipynb

Lastly, `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb` should appear in the
**time stamped** folder  Once all the MCMC BEAST 2 chains have finished running the copy of `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb` needs to be done manually. 
The progress of the MCMC BEAST 2 chains can be checked via openning the files `run-1.out`,
`run-2.out`, `run-3.out` and `run-4.out` or the command:
```bash
   squeue --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R" --me
```
The slurm jobs running MCMC BEAST 2 chains will be named `beast_full_sample_run_1`,
`beast_full_sample_run_2`, `beast_full_sample_run_3` and `beast_full_sample_run_4`.
Competing the manual run of the copy of `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb`
produces the output report notebook `BEAST_pype-Report.ipynb` and html version `BEAST_pype-Report.html`.

## Further Notes: 
* The slurm command `scancel -u YOUR-USERNAME` will cancel any slurm jobs you have 
   running, including the pipeline.