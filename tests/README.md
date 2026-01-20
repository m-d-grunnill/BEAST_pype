# Unit Tests

Unit testing is done using [PyTest](https://docs.pytest.org/en/stable/). To run all unit tests follow the instructions below. 
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
**Notes:** 
* NUMBER_OF_CPUS_TO_USE can be `logical` to run on all CPUs detected by the OS. BEWARE if working on a High Performance Cluster the detection method when using logical may detect more CPUs than allocated to you.
* It is not recommended that you run these tests from a HPC (like the one at NML). Doing so will cause the tests to fail with `RuntimeError("Kernel didn't respond in 60 seconds")`.