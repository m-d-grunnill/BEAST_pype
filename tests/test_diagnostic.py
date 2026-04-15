from click.testing import CliRunner
from beast_pype.cli import beast_pype
import os
import shutil
from pathlib import Path
from base_test_classes import _check_file_was_generated


def test_diagnostic(subtests):
    #copy needed files to a new directory to run the test
    src = "example_beast2_outputs_to_diagnose"
    dst = "Test_diagnostic"
    diagnostic_notebook = 'Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb'
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, dirs_exist_ok=True)
        else:
            # Ensure the destination directory exists
            os.makedirs(os.path.dirname(d), exist_ok=True)
            shutil.copy2(s, d)
    # run the tests.
    runner = CliRunner()
    run_args = ['diagnose-results', 'Generic', dst]
    result = runner.invoke(beast_pype, args=run_args)
    
    test_ran_ok_list = []
    with subtests.test(f"Check for error generation: {result.exc_info}"):
        test_ran_ok_list.append(result.exit_code == 0)
        assert test_ran_ok_list[-1]
    _check_file_was_generated(
        subtests=subtests,
        filename=diagnostic_notebook,
        directory=dst,
        test_ran_list=test_ran_ok_list
    )
    if all(test_ran_ok_list):
        shutil.rmtree(dst)