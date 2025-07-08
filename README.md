# BEAST_pype: A pipeline for phylodynamics using BEAST2.

## IMPORTANT NOTES
* With this early access version of beast_pype the beast_pype python package is
  installed with the `-e` flag. This means it is intsalled in an editable development mode.
  **MEANING IF YOU MAKE ANY CHANGES TO THE FILES IN  `src/beast_pype` THIS WILL
  AFFECT YOUR beast_pype INSTALLATION**
* Currently, beast_pype is reliant on [slurm](https://slurm.schedmd.com/overview.html). 
If you are running a beast_pype workflow on a machine without [slurm](https://slurm.schedmd.com/overview.html) it can be installed  via [conda](https://anaconda.org/conda-forge/slurm). However, the 
initial configuration of slurm may be quite involved.
* For ease of distribution reasons beast_pype uses the version of BEAST 2 that is available via conda, specifically [bioconda](https://anaconda.org/bioconda/beast2), 2.6 as 2025-Jun-25. Template beast2 xmls from other versions of may not work. BEAST 2.7.7 is available on the conda channel [millerjeremya](https://anaconda.org/millerjeremya/beast2). However, I tested this on a Linux OS (2025-06-04) and could not get the command line arguments to work.

## Initial Setup:
You should only need to carry out the following steps the first time you run the pipeline 
1. Clone the beast_pype repo: 
  ```bash
    git clone https://github.com/m-d-grunnill/BEAST_pype.git
 ```
OR:
  ```bash
    git clone git@github.com:m-d-grunnill/BEAST_pype.git
 ```

2. Install the beast_pype conda environment via conda or mamba:
```bash
conda env create -f requirements.yml
```
OR if you have access mamba (much faster):
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
wish to alter or add to (see  [slurm's sbatch documentation](https://slurm.schedmd.com/sbatch.html)) depending on a particular run of a BEAST_pype workflow.

**BEFORE RUNNING THE EXAMPLES BELOW** from the administrator of the HPC you
are using you will need to learn the name of the partition to use for your slurm jobs. 
Once you have acquired this you will need to alter line 11 of `run-workflow.slurm`  to read:
```bash
#SBATCH -p NAME_OF_PARTITION_TO_USE
```
Likewise, the parameter yml files below will need to 
have the line `partition: NAME_OF_PARTITION_TO_USE` added: 
* `parameters/Test-BDSKY-serial_full.yml`
* `parameters/Test-BDSKY-serial_no-initial-tree.yml`
* `parameters/Test-BDSKY-serial_xml-ready-to-go.yml`
* `parameters/Test-Generic_full.yml`
* `parameters/Test-Generic_no-initial-tree.yml`
* `parameters/Test-Generic_xml-ready-to-go.yml`

**NOTE** the values for `chain_length`,
`store_state_every`, `trace_log_every`, `tree_log_every`, `screen_log_every` and `store_state_every` in the files above are extremely low. As is the case for the xmls `example_data/COVID-19_BA.2.86/BDSky_serial_ready_to_go.xml` and `example_data/COVID-19_BA.2.86/Coalescent_Exponential_ready_to_go.xml`.  These low values are intended for testing that beast_pype workflows works on your set-up (taking a few minutes to run) and to provide an example of BEAST_pype's operation. **IT IS HIGHLY RECOMMEDED THAT** you up the values of `chain_length` (>=10,000,000),
`store_state_every` (>=1000), `trace_log_every` (>=1000), `tree_log_every` (>=1000), `screen_log_every` (>=1000) and `store_state_every`(>=1000) when properly running these workflows. 

### Running the BDSKY-serial Workflow



The [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb) along with other workflows
can be found in the `workflows` directory.  This workflow is intended for use with the Birth-Death Skyline (BDSKY)-serial model (see https://www.beast2.org/files/bdskytutorialv2.1.1.pdf). This will only work with a template
BDSky-serial xml or providing a ready-to-go BDSky-serial xml to use. When using a template
BDSky-serial xml, that template has the sequences from a fasta file inserted, (and depending on other
settings an initial tree added, the Rt dimensions and origin prior altered) before being run by BEAST 2.  A ready-to-go BDSky-serial xml will simply be run by BEAST 2. See parameters section at the top of [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb) for more detailed information on workflow parameters.  

The `parameters` folder contains example yml files that are used to run the different variants of the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) demonstrated below. You can alter the values in these files alter a copy of them for your own use.

#### Full Workflow

To run the full version o [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) from commandline via slurms `sbatch` enter the following command
```bash
cd LOCATION_YOU_CLONED_THE_BEAST_PYPE_REPO_TO
sbatch run-workflow.slurm workflows/BDSKY-serial.ipynb parameters/Test-BDSKY-serial_full.yml
```
This will run the workflow with the parameters set in `parameters/Test-BDSKY-serial_full.yml`. A description of the parameters used in [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb)
is at the top of the notebook. 

Running the above command will create the directory `example_runs_of_BSDKY`
containing a **time stamped** folder (format 'YYYY-MM-DD_hour-min-sec'). Shortly after the BDSKY-serial workflow is launched
a running copy of [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb) should appear.
The top code cell of this notebook will have a new code cell inserted below which replaces the default parameter values
in [BDSKY-serial workflow notebook](workflows/BDSKY-serial.ipynb) with the  values specified in
`parameters/Test-BDSKY-serial_full.yml`.  The 'time stamped' folder will gradually fill with running
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

<a id="BDSKY_skip_building"></a>
#### Skip building an initial tree and use BEAST 2's initial tree instead.

The full [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) creates an initial tree using IQTree and TreeTime (in phases 2i and 2ii, respectively) and uses the intial tree in generating a BEAST 2 xml from the template xml. This initial tree does speed up the convergence of MCMC chains when running BEAST 2. However, this is slightly against the spirit of MCMC analysis (see [BEAST 2 documentation](https://www.beast2.org/2014/07/28/all-about-starting-trees)).  [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) can skip building the initial tree if the line `use_initial_tree: False` is in the parameter yml (see `parameters/Test-BDSKY-serial_no-initial-tree.yml` for an example). An example of 
this can be launched via the command:

```bash
cd LOCATION_YOU_CLONED_THE_BEAST_PYPE_REPO_TO
sbatch run-workflow.slurm workflows/BDSKY-serial.ipynb parameters/Test-BDSKY-serial_no-initial-tree.yml
```

This variant of the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) will then run the in the same manner as the full version but miss Phases 2i and 2ii as such the Phase-2i-IQTree.ipynb,  
Phase-2ii-TreeTime-and-Down-Sampling.ipynb will not appear in the **time stamped**
folder (neither will any files associated with the running of IQtree or TreeTime).  

#### Using a ready to go xml.

You may have a BEAST BDSKY xml that you wish to run directly with
the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb).  This can be done by 
adding the variable `ready_to_go_xml` with the path to that BEAST BDSKY xml in
 the parameter yml file. The file `parameters/Test-BDSKY-serial_xml-ready-to-go.yml`
provides an example. To run this example use the command:

```bash
cd LOCATION_YOU_CLONED_THE_BEAST_PYPE_REPO_TO
sbatch run-workflow.slurm workflows/BDSKY-serial.ipynb parameters/Test-BDSKY-serial_xml-ready-to-go.yml
```

This variant of the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) will then run the in will miss phases 2-3
and copy the xml you provided into the **time stamped** folder and use it for running
the rest of the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb).

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
sbatch run-workflow.slurm workflows/Generic.ipynb parameters/Test-Generic_full.yml
```
This will run the workflow with the parameters set in `parameters/Test-Generic_full.yml`. 
You can alter them for usage or alter a copy of `parameters/Test-Generic_full.yml`.
A description of the parameters used in [Generic workflow notebook](workflows/Generic.ipynb)
is at the top of the notebook. 

Running the above command will create the directory `example_runs_of_BSDKY`
containing a **time stamped** folder (format 'YYYY-MM-DD_hour-min-sec'). Shortly after the Generic workflow is launched
a running copy of [Generic workflow notebook](workflows/Generic.ipynb) should appear.
The top code cell of this notebook will have a new code cell inserted below which replaces the default parameter values
in [Generic workflow notebook](workflows/Generic.ipynb) with the  values specified in
`parameters/Test-Generic_full.yml`.  The 'time stamped' folder will gradually fill with running
copies of the notebooks (in the order):
1. Phase-2i-IQTree.ipynb
2. Phase-2ii-TreeTime-and-Down-Sampling.ipynb
3. Phase-3-Gen-Generic-xml.ipynb
4. Phase-4-Running-BEAST.ipynb

Lastly, `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb` should appear in the
**time stamped** folder  Once all the MCMC BEAST 2 chains have finished running the copy of `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb` needs to be done manually. 
The progress of the MCMC BEAST 2 chains can be checked via opening the files `run-1.out`,
`run-2.out`, `run-3.out` and `run-4.out` or the command:
```bash
   squeue --format="%.18i %.9P %.30j %.8u %.8T %.10M %.9l %.6D %R" --me
```
The slurm jobs running MCMC BEAST 2 chains will be named `beast_full_sample_run_1`,
`beast_full_sample_run_2`, `beast_full_sample_run_3` and `beast_full_sample_run_4`.
Competing the manual run of the copy of `Phase-5-Diagnosing_Outputs_and_Generate_Report.ipynb`
produces the output report notebook `BEAST_pype-Report.ipynb` and html version `BEAST_pype-Report.html`.


#### Skip building an initial tree and use BEAST 2's initial tree instead.

As with the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) you can skip building an initial tree in the
[Generic workflow](workflows/Generic.ipynb) by adding `use_initial_tree: False` to the parameter yml
(see [above](#BDSKY_skip_building)). The `parameters/Test-Generic_no-initial-tree.yml` 
provides an example. To run this example use the command:

```bash
cd LOCATION_YOU_CLONED_THE_BEAST_PYPE_REPO_TO
sbatch run-workflow.slurm workflows/Generic.ipynb parameters/Test-Generic_no-initial-tree.yml
```

#### Using a ready to go xml.

The [Generic workflow](workflows/Generic.ipynb) like the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) can be run using a
beast 2 xml that is ready to go (does not need to be modified).  This can be done by 
adding the variable `ready_to_go_xml` with the path to that BEAST 2 xml in
 the parameter yml file. The file `parameters/Test-Generic_xml-ready-to-go.yml`
provides an example. To run this example use the command:

```bash
cd LOCATION_YOU_CLONED_THE_BEAST_PYPE_REPO_TO
sbatch run-workflow.slurm workflows/Generic.ipynb parameters/Test-Generic_xml-ready-to-go.yml
```

This variant of the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb) will then run the in will miss phases 2-3
and copy the xml you provided into the **time stamped** folder and use it for running
the rest of the [BDSKY-serial workflow](workflows/BDSKY-serial.ipynb).

## Further Notes: 
* The slurm command `scancel -u YOUR-USERNAME` will cancel any slurm jobs you have 
   running, including the pipeline.