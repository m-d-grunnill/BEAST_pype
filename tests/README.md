# Unit Tests

Unit testing is done using [PyTest](https://docs.pytest.org/en/stable/). To run all unit
tests:
1. If you have not installed the version of beast_pype you wish to test into the beast_pype environment do so via:
```bash
conda activate beast_pype
pip install .
```
2. Run the tests via the command below from the beast_pype conda environment (`conda activate beast_pype` to get into it).
```bash
pytest tests/
```