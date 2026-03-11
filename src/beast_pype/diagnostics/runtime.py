from subprocess import run
import io
import os
import pandas as pd


def _extracting_value_from_string(string, val_name):
    string = string.replace(val_name, '').strip()

    def is_part_of_number(character):
        if character.isdigit() or character == '.':
            return True
        else:
            return False

    return float(''.join(filter(is_part_of_number, string)))


def _extracting_values_from_txt_file(file_path, val_names):
    if isinstance(val_names, str):
        val_names = [val_names]
    file = open(file_path, "r")
    data = {}

    for line in file:
        for val_name in val_names:
            if line.startswith(val_name):
                if val_name in data:
                    raise AssertionError(val_name + ' is listed more than once in file.')
                data[val_name] = _extracting_value_from_string(line, val_name)

    return data


def get_beast_runtimes(directory,
                       outfile_startswith='run-with-seed-',
                       outfile_endswith='.out'):
    """
    Extract beast runtimes from out files.

    Parameters
    ----------
    directory: str
        Path to the directory containing the out files.
    outfile_startswith: str
        Out file ending suffix.
    outfile_endswith
        Out file prefix.

    Returns
    ---------
    pandas.DataFrame

    """
    time_records = []
    for file in os.listdir(directory):
        if file.startswith(outfile_startswith) and file.endswith(outfile_endswith):
            entry = {
                'BEAST_run': file.replace('.out', ''),
                **_extracting_values_from_txt_file(os.path.join(directory, file), 'Total calculation time: ')
            }
            time_records.append(entry)

    runtimes = pd.DataFrame(time_records)
    runtimes['run_time_D_H_M_S'] = pd.to_timedelta(runtimes['Total calculation time: '], unit='s')
    runtimes.rename(columns={'Total calculation time: ': 'run_time_seconds'}, inplace=True)
    return runtimes


def get_slurm_job_stats(job_ids):
    """
    Get statistics on Slurm jobs.

    Parameters
    ----------
    job_ids: list of strs
        List of job ids.

    Returns
    -------
    job_stats: pd.DataFrame
        DataFrame of slurm job statistics.

    """
    job_ids_request = ','.join([f"{entry}.batch" for entry in job_ids])
    request = f"sacct --jobs={job_ids_request} --format=JobID,AllocTres,Elapsed,CPUTime,TotalCPU,MaxRSS -p --delimiter='/t'"
    results = run(request, shell=True, capture_output=True, text=True)
    run_info_batch = pd.read_csv(io.StringIO(results.stdout), sep='/t')
    run_info_batch = run_info_batch.loc[:, ~run_info_batch.columns.str.startswith('Unnamed')]
    run_info_batch['JobID'] = run_info_batch['JobID'].str.replace('.batch', '')
    suffix_to_scientific = {'K': 'e3', 'M': 'e6', 'G': 'e9', 'T': 'e12'}
    maxrss = run_info_batch['MaxRSS']
    for suffix, scientific in suffix_to_scientific.items():
        maxrss = maxrss.str.replace(suffix, scientific)
    run_info_batch['Max RAM Used (GB)'] = maxrss.astype(float) / 1e9
    run_info_batch[['Allocated CPUs', 'Allocated RAM (GB)', 'Allocated Nodes']] = run_info_batch['AllocTRES'].str.split(
        ',', expand=True)
    run_info_batch['Allocated CPUs'] = run_info_batch['Allocated CPUs'].str.replace('cpu=', '').astype(int)
    run_info_batch['Allocated Nodes'] = run_info_batch['Allocated Nodes'].str.replace('node=', '').astype(int)
    run_info_batch['Allocated RAM (GB)'] = run_info_batch['Allocated RAM (GB)'].str.replace('mem=', '')
    run_info_batch['Allocated RAM (GB)'] = run_info_batch['Allocated RAM (GB)'].str.replace('G', '').astype(float)
    run_info_batch['Elapsed'] = pd.to_timedelta(run_info_batch['Elapsed'].str.replace('-', ' days '))
    run_info_batch['CPUTime'] = pd.to_timedelta(run_info_batch['CPUTime'].str.replace('-', ' days '))
    for i, val in enumerate(run_info_batch['TotalCPU'].str.count(':') == 1):
        if val:
            run_info_batch.loc[i, 'TotalCPU'] = '00:' + run_info_batch.loc[i, 'TotalCPU']
    if any(run_info_batch['TotalCPU'].str.contains('-', regex=False)):
        run_info_batch['TotalCPU'].str.replace('-', ' days ')
    run_info_batch['TotalCPU'] = pd.to_timedelta(run_info_batch['TotalCPU'])
    run_info_batch['CPU Efficiency (%)'] = 100 * run_info_batch['TotalCPU'] / run_info_batch['CPUTime']
    run_info_batch['RAM Efficiency (%)'] = 100 * run_info_batch['Max RAM Used (GB)'] / run_info_batch[
        'Allocated RAM (GB)']

    job_ids_request = ','.join(job_ids)
    request = f"sacct --jobs={job_ids_request} --format=JobID,JobName,Timelimit -p --delimiter='/t'"
    results = run(request, shell=True, capture_output=True, text=True)
    run_info = pd.read_csv(io.StringIO(results.stdout), sep='/t')
    run_info = run_info[run_info['JobID'].isin(job_ids)]
    run_info = run_info.loc[:, ~run_info.columns.str.startswith('Unnamed')]
    run_info['Timelimit'] = pd.to_timedelta(run_info['Timelimit'].str.replace('-', ' days '))

    job_stats = run_info.merge(run_info_batch, on='JobID')
    job_stats['Timelimit Used %'] = 100 * job_stats['Elapsed'] / job_stats['Timelimit']
    job_stats = job_stats[['JobID', 'JobName', 'Elapsed', 'Timelimit', 'Timelimit Used %', 'Allocated Nodes',
                           'Allocated CPUs', 'TotalCPU', 'CPUTime', 'CPU Efficiency (%)', 'Max RAM Used (GB)',
                           'Allocated RAM (GB)', 'RAM Efficiency (%)']]
    return job_stats