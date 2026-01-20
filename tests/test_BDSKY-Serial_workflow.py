from base_test_classes import SimpleWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-BDSKY-Serial-xml.ipynb'
workflow = 'BDSKY-Serial'

class TestBDSKYSerialFull(SimpleWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_full.yml'
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialChangeTimes(SimpleWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_change_times.yml'
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialNoInitialTree(SimpleWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_no-initial-tree.yml'
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialXmlReadyToGo(SimpleWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/BDSKY-Serial_xml-ready-to-go.yml'
    workflow = workflow
    variation = 'xml ready-to-go'
    xml_generation_notebook = xml_generation_notebook