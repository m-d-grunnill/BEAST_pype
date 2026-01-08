from base_test_classes import ComparativeWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-Generic-xml.ipynb'
workflow = 'Generic-Comparative'

class TestGenericComparativeFull(ComparativeWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/Generic-Comparative_full.yml'
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

class TestGenericComparativeNoInitialTree(ComparativeWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/Generic-Comparative_no-initial-tree.yml'
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook