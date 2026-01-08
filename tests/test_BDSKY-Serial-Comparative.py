from base_test_classes import ComparativeWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-BDSKY-Serial-xml.ipynb'
workflow = 'BDSKY-Serial-Comparative'

class TestBDSKYSerialComparativeFull(ComparativeWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial-Comparative_full.yml'
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialComparativeChangeTimes(ComparativeWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial-Comparative_change_times.yml'
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialComparativeNoInitialTree(ComparativeWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial-Comparative_no-initial-tree.yml'
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook