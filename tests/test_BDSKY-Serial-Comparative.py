from base_test_classes import ComparativeWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-BDSKY-Serial-xml.ipynb'
workflow = 'BDSKY-Serial-Comparative'

class TestBDSKYSerialComparativeFull(ComparativeWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_Comparative_Full'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial-Comparative_full.yml'
    workflow = workflow
    initial_tree_building =True
    downsampling = True
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialComparativeNoDownsampling(ComparativeWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_Comparative_No_Downsampling'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial-Comparative_no_downsampling.yml'
    workflow = workflow
    initial_tree_building =True
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook


class TestBDSKYSerialComparativeChangeTimes(ComparativeWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_Comparative_Change_Times'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial-Comparative_partitions.yml'
    workflow = workflow
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialComparativeNoInitialTree(ComparativeWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_Comparative_No_Initial_Tree'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial-Comparative_no-initial-tree.yml'
    workflow = workflow
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook