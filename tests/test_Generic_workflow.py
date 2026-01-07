from base_test_classes import SimpleWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-Generic-xml.ipynb'
workflow = 'Generic'

class TestGenericFull(SimpleWorkflowVariationTest):
    parameters = {
        'collection_date_field': 'sample collection date',
        'fasta_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_alignment.fasta',
        'initial_tree_type': 'Temporal',
        'down_sample_to': 100,
        'log_file_basename': 'BEAST',
        'metadata_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_metadata.csv',
        'number_of_beast_runs': 3,
        'overall_save_dir': 'Local-Test-Generic_full',
        'sample_id_field': 'specimen collector sample ID',
        'template_xml_path': 'template_beast_xmls/Coalescent_Exponential_COVID-19_template.xml',
        'chain_length': 5000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000
    }
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

class TestGenericNoInitialTree(SimpleWorkflowVariationTest):
    parameters = {
        'collection_date_field': 'sample collection date',
        'fasta_path': 'example_data/COVID-19_JN.1/down_sampled_sequences.fasta',
        'use_initial_tree': False,
        'log_file_basename': 'BEAST',
        'metadata_path': 'example_data/COVID-19_JN.1/down_sampled_metadata.csv',
        'number_of_beast_runs': 3,
        'overall_save_dir': 'Local-Test-Generic_no-initial-tree',
        'sample_id_field': 'specimen collector sample ID',
        'template_xml_path': 'template_beast_xmls/Coalescent_Exponential_COVID-19_template.xml',
        'chain_length': 5000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000
    }
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestGenericXmlReadyToGo(SimpleWorkflowVariationTest):
    parameters = {
        'number_of_beast_runs': 3,
        'overall_save_dir': 'Local-Test-Generic_xml-ready-to-go',
        'ready_to_go_xml': 'example_data/COVID-19_JN.1/Coalescent_Exponential_ready_to_go.xml'
    }
    workflow = workflow
    variation = 'xml ready-to-go'
    xml_generation_notebook = xml_generation_notebook

