from base_test_classes import ComparativeWorkflowVariationTest

xml_generation_notebook = 'Phase-3-Gen-BDSKY-Serial-xml.ipynb'
workflow = 'COVID-Strain-Surveillance'

class TestCOVIDStrainSurveillance(ComparativeWorkflowVariationTest):
    name = 'Test_COVID_Strain_Surveillance'
    parameters_path = 'parameters/locally_run_examples/COVID-Strain-Surveillance.yml'
    workflow = workflow
    variation = 'full'
    xml_generation_notebook = xml_generation_notebook

    @property
    def xml_set_labels(self):
        return [f'VOI_{strain}' for strain in self.parameters['voi_strains']] + [f"DR_{self.parameters['dr_strain']}"]

