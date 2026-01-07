from base_test_classes import SimpleWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-BDSKY-Serial-xml.ipynb'
workflow = 'BDSKY-Serial'

class TestBDSKYSerialFull(SimpleWorkflowVariationTest):
    parameters = {
        'collection_date_field': 'sample collection date',
        'fasta_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_alignment.fasta',
        'down_sample_to': 100,
        'initial_tree_type': 'Temporal',
        'log_file_basename': 'BDSKY_serial',
        'metadata_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_metadata.csv',
        'number_of_beast_runs': 3,
        'origin_start_addition': 0.02737850787132101,
        'origin_upper_addition': 2,
        'overall_save_dir': 'Local-Test-BDSKY-Serial_full',
        'rt_dims': 3,
        'zero_sampling_before_first_sample': True,
        'sample_id_field': 'specimen collector sample ID',
        'template_xml_path': 'template_beast_xmls/BDSKY_serial_COVID-19_template.xml',
        'chain_length': 50000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000
    }
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialChangeTimes(SimpleWorkflowVariationTest):
    parameters = {
        'collection_date_field': 'sample collection date',
        'fasta_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_alignment.fasta',
        'use_initial_tree': False,
        'log_file_basename': 'BDSKY_serial',
        'metadata_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_metadata.csv',
        'number_of_beast_runs': 3,
        'overall_save_dir': 'Local-Test-BDSKY_change_times',
        'rt_changes': {'unit': 'months', 'every': 1},
        'sampling_prop_changes': {'unit': 'months', 'every': 1},
        'zero_sampling_before_first_sample': True,
        'sample_id_field': 'specimen collector sample ID',
        'template_xml_path': 'template_beast_xmls/BDSKY_serial_COVID-19_template.xml',
        'chain_length': 50000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000
    }
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialNoInitialTree(SimpleWorkflowVariationTest):
    parameters =  {
        'collection_date_field': 'sample collection date',
        'fasta_path': 'example_data/COVID-19_JN.1/down_sampled_sequences.fasta',
        'use_initial_tree': False,
        'log_file_basename': 'BDSKY_serial',
        'metadata_path': 'example_data/COVID-19_JN.1/down_sampled_metadata.csv',
        'number_of_beast_runs': 3,
        'overall_save_dir': 'Local-Test-BDSKY-Serial_no-initial-tree',
        'rt_dims': 3,
        'zero_sampling_before_first_sample': True,
        'sample_id_field': 'specimen collector sample ID',
        'template_xml_path': 'template_beast_xmls/BDSKY_serial_COVID-19_template.xml',
        'chain_length': 50000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000
    }
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialXmlReadyToGo(SimpleWorkflowVariationTest):
    parameters = {
        'number_of_beast_runs': 3,
        'overall_save_dir': 'Local-Test-BDSKY-Serial_xml-ready-to-go',
        'ready_to_go_xml': 'example_data/COVID-19_JN.1/BDSky_test_ready_to_go.xml'
    }
    workflow = workflow
    variation = 'xml ready-to-go'
    xml_generation_notebook = xml_generation_notebook