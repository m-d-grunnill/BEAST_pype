# Regression Tests

Testing is done using [PyTest](https://docs.pytest.org/en/stable/). To run all tests follow the instructions below. 
Note these instructions assume you are at the root of the repository.
1. If you have not installed the version of beast_pype you wish to test into the beast_pype environment do so via:
```bash
conda activate beast_pype
pip install .
```


2. Run the tests via the command below from the beast_pype conda environment (`conda activate beast_pype` to get into it).
```bash
pytest -n NUMBER_OF_CPUS_TO_USE tests/
```
**Note** NUMBER_OF_CPUS_TO_USE can be `logical` to run on all CPUs detected by the OS. BEWARE if working on a High Performance Cluster the detection method when using logical may detect more CPUs than allocated to you.

# Running the tests against a specific Jupyter kernelspec

The workflow tests execute notebooks using a Jupyter kernel. By default they use the
`beast_pype` kernelspec. If you have installed beast_pype into a differently named
environment (and registered a kernelspec for it), you can point the tests at that
kernelspec with the `--kernel-name` option:
```bash
pytest -n NUMBER_OF_CPUS_TO_USE tests/ --kernel-name KERNELSPEC_NAME
```
For example, to run the tests with 31 CPUs against a development environment whose
kernelspec is named `dev_beast_pype`:
```bash
pytest -n 31 tests/ --kernel-name dev_beast_pype
```
**Note** `KERNELSPEC_NAME` must be the name of an installed kernelspec (not the conda
environment path). You can list the kernelspecs available on your machine with:
```bash
jupyter kernelspec list
```
If a kernelspec for your environment does not yet exist, you can register one from
within that environment via:
```bash
python -m ipykernel install --user --name KERNELSPEC_NAME --display-name "KERNELSPEC_NAME"
```

# A note on running these tests with little RAM 

If you do not have much RAM tests may well fail with the error `RuntimeError("Kernel didn't respond in 60 seconds")`. There are 2 solutions to this:
1. Run the tests serially via the command `pytest -n 1 tests/`.
2. Within the file `base_test_classes.py` alter the class `WorkflowVariationTest's` `start_timeout` attribute to be 600 (i.e. 10 minutes) or higher.

Unfortunately, these solutions can greatly increase the run time of all these tests.
