from base_test_classes import SimpleWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-Generic-xml.ipynb'
workflow = 'Generic'

class TestGenericFull(SimpleWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/Generic_full.yml'
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

class TestGenericNoInitialTree(SimpleWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/Generic_no-initial-tree.yml'
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestGenericXmlReadyToGo(SimpleWorkflowVariationTest):
    parameters_path = 'parameters/locally_run_examples/Generic_xml-ready-to-go.yml'
    workflow = workflow
    variation = 'xml ready-to-go'
    xml_generation_notebook = xml_generation_notebook