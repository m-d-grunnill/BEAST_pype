from base_test_classes import SimpleWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-BDSKY-Serial-xml.ipynb'
workflow = 'BDSKY-Serial'

class TestBDSKYSerialFull(SimpleWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_Full'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_full.yml'
    workflow = workflow
    initial_tree_building = True
    downsampling = True
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialNoDownsampling(SimpleWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_No_Downsampling'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_no_downsampling.yml'
    workflow = workflow
    initial_tree_building = True
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook


class TestBDSKYSerialPartitions(SimpleWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_Partitions'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_partitions.yml'
    workflow = workflow
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialNoInitialTree(SimpleWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_No_Initial_Tree'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_no-initial-tree.yml'
    workflow = workflow
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialXmlReadyToGo(SimpleWorkflowVariationTest):
    name = 'Test_BDSKY-Serial_XML_Ready_To_Go'
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_xml-ready-to-go.yml'
    workflow = workflow
    xml_generation_notebook = xml_generation_notebook