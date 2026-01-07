from base_test_classes import ComparativeWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-BDSKY-Serial-xml.ipynb'
workflow = 'BDSKY-Serial-Comparative'

class TestBDSKYSerialComparativeFull(ComparativeWorkflowVariationTest):
    parameters = {'overall_save_dir': 'Local-Test-BDSKY-Serial-Comparative_full',
                  'template_xml_path': 'template_beast_xmls/BDSKY_serial_COVID-19_template.xml',
                  'fasta_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_alignment.fasta',
                  'metadata_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_metadata.csv',
                  'sample_id_field': 'specimen collector sample ID',
                  'collection_date_field': 'sample collection date',
                  'xml_set_definitions': {'JN.1': "`lineage name` == 'JN.1'",
                                          'JN.1.1': "`lineage name` == 'JN.1.1'"},
                  'use_initial_tree': True,
                  'down_sample_to': 100,
                  'rt_dims': 3,
                  'zero_sampling_before_first_sample': True,
                  'log_file_basename': 'BDSKY_Serial',
                  'chain_length': 50000,
                  'screen_log_every': 5000,
                  'store_state_every': 1000,
                  'trace_log_every': 1000,
                  'tree_log_every': 1000,
                  'number_of_beast_runs': 3}
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialComparativeChangeTimes(ComparativeWorkflowVariationTest):
    parameters = {
        'overall_save_dir': 'Local-Test-BDSKY-Serial-Comparative_change_times',
        'template_xml_path': 'template_beast_xmls/BDSKY_serial_COVID-19_template.xml',
        'fasta_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_alignment.fasta',
        'metadata_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_metadata.csv',
        'sample_id_field': 'specimen collector sample ID',
        'collection_date_field': 'sample collection date',
        'use_initial_tree': False,
        'xml_set_definitions': {'JN.1': "`lineage name` == 'JN.1'",
                                'JN.1.1': "`lineage name` == 'JN.1.1'"},
        'rt_changes': {'unit': 'months', 'every': 1},
        'sampling_prop_changes': {'unit': 'months', 'every': 1},
        'zero_sampling_before_first_sample': True,
        'log_file_basename': 'BDSKY_serial',
        'chain_length': 50000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000,
        'number_of_beast_runs': 3
    }
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook

class TestBDSKYSerialComparativeNoInitialTree(ComparativeWorkflowVariationTest):
    parameters = {
        'overall_save_dir': 'Local-Test-BDSKY-Serial-Comparative_no-initial-tree',
        'template_xml_path': 'template_beast_xmls/BDSKY_serial_COVID-19_template.xml',
        'fasta_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_alignment.fasta',
        'metadata_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_metadata.csv',
        'sample_id_field': 'specimen collector sample ID',
        'collection_date_field': 'sample collection date',
        'xml_set_definitions': {'JN.1': "`lineage name` == 'JN.1'",
                                'JN.1.1': "`lineage name` == 'JN.1.1'"},
        'use_initial_tree': False,
        'rt_dims': 3,
        'zero_sampling_before_first_sample': True,
        'log_file_basename': 'BDSKY_Serial',
        'chain_length': 50000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000,
        'number_of_beast_runs': 3
    }
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook
