from base_test_classes import ComparativeWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-Generic-xml.ipynb'
workflow = 'Generic-Comparative'

class TestGenericComparativeFull(ComparativeWorkflowVariationTest):
    parameters = {
        'collection_date_field': 'sample collection date',
        'fasta_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_alignment.fasta',
        'down_sample_to': 100,
        'metadata_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_metadata.csv',
        'number_of_beast_runs': 3,
        'overall_save_dir': 'Local-Test-Generic-Comparative_full',
        'log_file_basename': 'BEAST',
        'sample_id_field': 'specimen collector sample ID',
        'template_xml_path': 'template_beast_xmls/Coalescent_Exponential_COVID-19_template.xml',
        'chain_length': 5000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000,
        'xml_set_definitions': {'JN.1': "`lineage name` == 'JN.1'",
                                'JN.1.1': "`lineage name` == 'JN.1.1'"}
    }
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

class TestGenericComparativeNoInitialTree(ComparativeWorkflowVariationTest):
    parameters = {
        'collection_date_field': 'sample collection date',
        'fasta_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_alignment.fasta',
        'use_initial_tree': False,
        'metadata_path': 'example_data/COVID-19_JN.1/VirusSeq_JN1_metadata.csv',
        'number_of_beast_runs': 3,
        'overall_save_dir': 'Local-Test-Generic-Comparative_no-initial-tree',
        'log_file_basename': 'BEAST',
        'sample_id_field': 'specimen collector sample ID',
        'template_xml_path': 'template_beast_xmls/Coalescent_Exponential_COVID-19_template.xml',
        'chain_length': 5000,
        'screen_log_every': 5000,
        'store_state_every': 1000,
        'trace_log_every': 1000,
        'tree_log_every': 1000,
        'xml_set_definitions': {'JN.1': "`lineage name` == 'JN.1'",
                                'JN.1.1': "`lineage name` == 'JN.1.1'"}
    }
    workflow = workflow
    variation = 'no initial tree'
    xml_generation_notebook = xml_generation_notebook
