from base_test_classes import SimpleWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-Generic-xml.ipynb'
workflow = 'Generic'

class TestGenericFull(SimpleWorkflowVariationTest):
    name = 'Test_Generic_Full'
    parameters_path = 'parameters/locally_run_examples/Generic_full.yml'
    workflow = workflow
    initial_tree_building = True
    xml_generation = True
    downsampling = True
    xml_generation_notebook = xml_generation_notebook

class TestGenericNoDownsampling(SimpleWorkflowVariationTest):
    name = 'Test_Generic_No_Downsampling'
    parameters_path = 'parameters/locally_run_examples/Generic_no_downsampling.yml'
    workflow = workflow
    initial_tree_building = True
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestGenericNoInitialTree(SimpleWorkflowVariationTest):
    name = 'Test_Generic_No_Initial_Tree'
    parameters_path = 'parameters/locally_run_examples/Generic_no-initial-tree.yml'
    workflow = workflow
    xml_generation = True
    xml_generation_notebook = xml_generation_notebook

class TestGenericXmlReadyToGo(SimpleWorkflowVariationTest):
    name = 'Test_Generic_XML_Ready_To_Go'
    parameters_path = 'parameters/locally_run_examples/Generic_xml-ready-to-go.yml'
    workflow = workflow
    xml_generation_notebook = xml_generation_notebook