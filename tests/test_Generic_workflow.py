from base_test_classes import SimpleWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-Generic-xml.ipynb'
workflow = 'Generic'

class TestGenericInitialTree(SimpleWorkflowVariationTest):
    name = 'Test_Generic_Initial_Tree'
    parameters_path = 'parameters/locally_run_examples/Generic_InitialTree.yml'
    workflow = workflow
    variation = 'initial tree building'
    xml_generation_notebook = xml_generation_notebook

class TestGenericNoInitialTree(SimpleWorkflowVariationTest):
    name = 'Test_Generic_No_Initial_Tree'
    parameters_path = 'parameters/locally_run_examples/Generic_no-initial-tree.yml'
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestGenericXmlReadyToGo(SimpleWorkflowVariationTest):
    name = 'Test_Generic_XML_Ready_To_Go'
    parameters_path = 'parameters/locally_run_examples/Generic_xml-ready-to-go.yml'
    workflow = workflow
    variation = 'xml ready-to-go'
    xml_generation_notebook = xml_generation_notebook