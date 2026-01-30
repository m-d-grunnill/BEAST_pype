from base_test_classes import ComparativeWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-Generic-xml.ipynb'
workflow = 'Generic-Comparative'

class TestGenericComparativeFull(ComparativeWorkflowVariationTest):
    name = 'Test_Generic_Comparative_Full'
    parameters_path = 'parameters/locally_run_examples/Generic-Comparative_full.yml'
    workflow = workflow
    initial_tree_building = True
    downsampling = True
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestGenericComparativeNoDownsampling(ComparativeWorkflowVariationTest):
    name = 'Test_Generic_Comparative_No_Downsampling'
    parameters_path = 'parameters/locally_run_examples/Generic-Comparative_no_downsampling.yml'
    workflow = workflow
    initial_tree_building = True
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestGenericComparativeNoInitialTree(ComparativeWorkflowVariationTest):
    name = 'Test_Generic_Comparative_No_Initial_Tree'
    parameters_path = 'parameters/locally_run_examples/Generic-Comparative_no-initial-tree.yml'
    workflow = workflow
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook