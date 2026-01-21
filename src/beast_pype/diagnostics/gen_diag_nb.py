import nbformat as nbf
import importlib.resources as importlib_resources
import os
from papermill.iorw import load_notebook_node, write_ipynb
from papermill.parameterize import parameterize_notebook
from papermill.inspection import _infer_parameters
from papermill.iorw import load_notebook_node



workflow_modules = importlib_resources.path('beast_pype', 'workflow_modules')
reports_path = importlib_resources.path('beast_pype', 'report_templates')
available_reports = [file.replace('.ipynb', '') for file in os.listdir(reports_path) if file.endswith('.ipynb')]

def gen_xml_set_diag_notebook(save_dir,
                              directories_to_exclude=None):
    """
    Generates rest of Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb.

    Parameters
    ----------
    save_dir: str
        Path to save copy of Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb.
        Each subdirectory of save_dir will be assumed to contain data of a xml_set
         and used in generating the new Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb.
    directories_to_exclude:
        Any directories to exclude from save_dir when generating new
        Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb.

    Returns
    -------
    None
    """
    if directories_to_exclude is None:
        directories_to_exclude = ['.ipynb_checkpoints']
    directories = [directory for directory in os.listdir(save_dir)
                   if os.path.isdir(os.path.join(save_dir, directory)) and
                   directory not in directories_to_exclude]
    workflows_modules = importlib_resources.path('beast_pype', 'workflow_modules')
    diag_notebook_path = f'{workflows_modules}/Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb'
    diag_notebook = nbf.read(diag_notebook_path, as_version=4)
    for xml_set in directories:
        diag_notebook['cells'] += [
            nbf.v4.new_markdown_cell(
                f"## Diagnosing XML Set: {xml_set}\n" +
                "### Loading data"),
            nbf.v4.new_code_cell(
                f"if os.path.exists('{xml_set}/beast_outputs'):\n" +
                f"\tbeast_outputs_path = '{xml_set}/beast_outputs'\n" +
                "else:\n" +
                f"\tbeast_outputs_path = '{xml_set}'\n\n"
                "sample_diag = BEASTDiag(beast_outputs_path)\n"+
                f"beast_outputs_paths['{xml_set}'] = beast_outputs_path"
            ),
            nbf.v4.new_markdown_cell(
                "## Selecting burn-in and Chains to Remove\n\n" +
                "Activating the cell below will generate an interactive widget. Widgets parts:\n" +
                "* Top interactive part: this allows you to select for a different burn-in and remove chains and select the parameters used in the rest of the widget.,\n" +
                "* Middle display: KDE and trace plots, see [arviz.plot_trace documentation](https://python.arviz.org/en/stable/api/generated/arviz.plot_trace.html#arviz.plot_trace).\n" +
                "* Bottom display: A table of statistics regarding the traces, see [arviz.summary documentation](https://python.arviz.org/en/stable/api/generated/arviz.summary.html#arviz.summary). Regarding these statistics:\n" +
                "\t* Ideally the ESSs should be >= 200, see [arviz.ess documentation](https://python.arviz.org/en/stable/api/generated/arviz.ess.html#arviz.ess).\n" +
                "\t* Ideally the r_hat should be close fo 1, see [arviz.rhat documentation](https://python.arviz.org/en/stable/api/generated/arviz.rhat.html#arviz.rhat).\n" +
                "\t* Markov Chain Standard Error MCSEs, see [arviz.mcse](https://python.arviz.org/en/stable/api/generated/arviz.mcse.html#arviz.mcse).\n\n"+
                "After making your selection click on the cell below the widget and then keep pressing shift+enter to carry on with the rest of the cells in this notebook."
            ),
            nbf.v4.new_code_cell(
                "sample_diag_widget = sample_diag.generate_widget(parameters_displayed=4)\n" +
                "sample_diag_widget"),
            nbf.v4.new_code_cell(
                f"pipeline_run_info['Chains Used']['{xml_set}'] = deepcopy(sample_diag.selected_chains)\n" +
                f"pipeline_run_info['Burn-In']['{xml_set}'] = deepcopy(sample_diag.burinin_percentage)\n" +
                f"phase_5i_params = sample_diag.merging_outputs_params(output_path=outputs_and_reports_dir, xml_set='{xml_set}')\n" +
                "phase_5i_log = execute_notebook(input_path=f'{workflow_modules}/Phase-5i-Merge-BEAST-outputs.ipynb',\n" +
                 f"output_path=save_dir + '/{xml_set}/Phase-5i-Merge-BEAST-outputs.ipynb',\n" +
                "\t\t\t\t\t\t\t\tparameters=phase_5i_params,\n" +
                "\t\t\t\t\t\t\t\tprogress_bar=True,\n" +
                "\t\t\t\t\t\t\t\tnest_asyncio=True)"
            )
        ]

    diag_notebook['cells'] += [
        nbf.v4.new_markdown_cell(
            "## Update the pipeline_run_info yaml."),
        nbf.v4.new_code_cell(
            "with open(f'{save_dir}/pipeline_run_info.yml', 'w') as fp:\n\tyaml.dump(pipeline_run_info, fp, sort_keys=True)\nfp.close()"
        ),
        nbf.v4.new_markdown_cell('## Get BEAST 2 runtimes'),
        nbf.v4.new_code_cell(
            "runtimes_dfs = []\n"+
            "job_stats_dfs = []\n"
            "for xml_set, path in beast_outputs_paths.items():\n" +
            "\tentry_df = get_beast_runtimes(path, outfile_startswith='run-with-seed-', outfile_endswith='.out')\n"
            "\tentry_df['xml_set'] = xml_set\n" +
            "\truntimes_dfs.append(entry_df)\n" +
            "\tif os.path.isfile(f'{path}/slurm_job_ids.txt'):\n" +
            "\t\tjobs_df = pd.read_csv(f'{path}/slurm_job_ids.txt', sep=';')\n" +
            "\t\tjobs_df['JobID']=jobs_df['JobID'].astype(str)\n" +
            "\t\tstats_df = get_slurm_job_stats(jobs_df['JobID'].to_list())\n" +
            "\t\tjob_stats_entry = jobs_df.merge(stats_df, on='JobID')\n" +
            "\t\tjob_stats_entry['xml_set'] = xml_set\n" +
            "\t\tjob_stats_dfs.append(job_stats_entry)\n" +
            "runtimes_df = pd.concat(runtimes_dfs)\n"+
            "runtimes_df.to_csv(f'{outputs_and_reports_dir}/BEAST_runtimes.csv', index=False)\n"+
            "if job_stats_dfs:\n"+
            "\tjob_stats_df = pd.concat(job_stats_dfs)\n" +
            "\tjob_stats_df.to_csv(f'{outputs_and_reports_dir}/BEAST_slurm_stats.csv', index=False)"
        ),
        nbf.v4.new_markdown_cell(
            '## Generate output Report\n' +
            'Now you can now move on to visualising outputs from BEAST using a report template.'
        ),
        nbf.v4.new_code_cell(
            "report_params = {'save_dir': outputs_and_reports_dir, 'beast_xml_path':beast_xml_path, 'xml_set_label': xml_set_label}\n" +
            "output_report_path = f'{outputs_and_reports_dir}/BEAST_pype-Report.ipynb'\n" +
            "add_unreported_outputs(report_template, outputs_and_reports_dir, output_report_path, xml_set_comparisons=True)\n" +
            "output = execute_notebook(\n" +
            "\tinput_path=output_report_path,\n" +
            "\toutput_path=output_report_path,\n" +
            "\tparameters=report_params,\n" +
            "\tkernel_name=kernel_name,\n" +
            "\tprogress_bar=True)"),
        nbf.v4.new_markdown_cell(
            "### Convert Output Report from Jupyter Notebook to Notebook\n" +
            "This also removes code cells."
        ),
        nbf.v4.new_code_cell(
            "%%bash -l -s {output_report_path}\n" +
            "source activate beast_pype\n" +
            "jupyter nbconvert --to html --no-input $@"
        ),
        nbf.v4.new_markdown_cell(
            "## Produce MCC trees\n\n"+
            "This can take sometime. So you can look at the report whilst waiting this is done last."
        ),
        nbf.v4.new_code_cell(
            "gen_mcc_notebook(outputs_and_reports_dir, 'Phase-5ii-Gen-MCC-Trees.ipynb')\n" +
            "mcc_tree_output = execute_notebook(input_path='Phase-5ii-Gen-MCC-Trees.ipynb',\n" +
            "output_path='Phase-5ii-Gen-MCC-Trees.ipynb',\n" +
            "progress_bar=True)\n"
        )
    ]

    with open(f'{save_dir}/Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb', 'w') as f:
        nbf.write(diag_notebook, f)

def gen_beast_diagnostic_nb(beast_outputs,
                            report_template,
                            kernel_name='beast_pype',
                            beast_xml_path=None, **kwargs):
    """
    Generate notebook for Diagnosing BEAST 2 MCMC runs & generating a report.

    Parameters
    ----------
    beast_outputs: str
        A valid path containing output `.log` and `.tree` files from running BEAST 2. The
        Phase 5 diagnostic report is generated within this directory.
        On how BEAST 2 outputs should be organised:
        *  If report_template is for simple reports ('BDSKY-Serial',
            or 'Generic') this should contain BEAST 2 outputs you wish to diagnose.
        *   If report_template is for comparative reports ('BDSKY-Serial-Comparative',
            'COVID-Strain-Surveillance' or 'Generic-Comparative') this should contain
            subdirectories each containing the BEAST 2 outputs you wish to diagnose and compare.
        Note in both cases beast outputs can be contained in a subdirectory labeled `beast_outputs`.
    report_template: str
        Name of a valid report template to use to generate report.
    kernel_name: str, default 'beast_pype'
        Name of Jupyter python kernel_name to use when running diagnostic & report template notebooks.
    beast_xml_path: str or dict of strings, optional
        Path(s) to BEAST 2 xml used to generate BEAST outputs. If not given, the
        location described in beast_outputs is checked for the file `beast.xml`
        and that is used instead.
        If report_template is for comparative reports ('BDSKY-Serial-Comparative',
        'COVID-Strain-Surveillance' or 'Generic-Comparative') this should be a dictionary
        with keys being the names of xml_sets (names of subdirectories described in
        beast_outputs) values being the paths to the corresponding BEAST 2 xml file.



    """
    if report_template not in available_reports:
        raise ValueError('The value for report_template must be one of the following: \n{}'.format(', '.join(available_reports)))
    report_template_path = f"{reports_path}/{report_template}.ipynb"
    report_nb = load_notebook_node(report_template_path)
    accepted_parameters = set([parameter.name
                           for parameter in _infer_parameters(report_nb)])
    invalid_params = set(kwargs.keys()) - accepted_parameters
    if invalid_params:
        raise ValueError(
            f'The following parameters that have been supplied for {report_template} as kwargs are invalid: \n {", ".join(invalid_params)}')
    parameters = kwargs
    comparative_reports = ['BDSKY-Serial-Comparative',
                           'COVID-Strain-Surveillance',
                           'Generic-Comparative']
    if report_template in comparative_reports:
        if beast_xml_path is None:
            beast_xml_path = {}
            for entry in os.listdir(beast_outputs):
                if os.path.isdir(f"{beast_outputs}/{entry}") and entry != '.ipynb_checkpoints':
                    beast_xml_path_value = f"{entry}/beast.xml"
                    if not os.path.isfile(f"{beast_outputs}/{beast_xml_path_value}"):
                        beast_xml_path_value = f"{entry}/beast_outputs/beast.xml"
                    beast_xml_path[entry] = beast_xml_path_value
        gen_xml_set_diag_notebook(beast_outputs, directories_to_exclude=['.ipynb_checkpoints'])
        diagnostic_save_name = f"{beast_outputs}/Phase-5-Diagnosing-XML-sets-and-Generate-Report.ipynb"
        diagnostic_nb = load_notebook_node(diagnostic_save_name)
    else:
        if beast_xml_path is None:
            beast_xml_path = "beast.xml"
            if not os.path.isfile(f"{beast_outputs}/{beast_xml_path}"):
                beast_xml_path = "beast_outputs/beast.xml"
        diagnostic_save_name = f'{beast_outputs}/Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb'
        diagnostic_nb = load_notebook_node(f'{workflow_modules}/Phase-5-Diagnosing-Outputs-and-Generate-Report.ipynb')
    parameters['beast_xml_path'] = beast_xml_path
    parameters['report_template'] = report_template_path
    parameters['kernel_name'] = kernel_name
    diagnostic_nb['metadata']['kernelspec']['name'] = kernel_name
    diagnostic_nb['metadata']['kernelspec']['display_name'] = kernel_name
    diagnostic_nb = parameterize_notebook(diagnostic_nb, parameters, kernel_name)
    write_ipynb(diagnostic_nb, diagnostic_save_name)

