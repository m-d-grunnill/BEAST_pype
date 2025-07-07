import time
import os
from subprocess import run
import io
import pandas as pd


def check_file_for_phrase(file_path,
                          phrase='slurm_job_complete',
                          wait_time=60):
    """
    If file exists check file for presence of phrase. Not found wait_time in seconds.

    Used to check a slurm .out file if the slurm job is complete,
    provided slurm sbatch ends with `; echo slurm_job_complete`.

    Parameters
    ----------
    file_path: str
       Path to text file (slurm .out file).
    phrase: str (default='slurm_job_complete')
       Phrase to check for file present at file_path.
    wait_time : int (default=60)
       Time in seconds to wait between each check.

    Returns
    -------
    None

    """
    # Section is silenced as the slurm job is not started imediatly.
    #if not os.path.exists(file_path):
    #    raise FileNotFoundError(f'file_path "{file_path}" does not exist.')
    #if not os.path.isfile(file_path):
    #    raise FileNotFoundError(f'file_path "{file_path}" is not a file.')
        
    complete_phrase_found = False
    while not complete_phrase_found:
        if os.path.isfile(file_path):
            with open(file_path) as file:
                s = file.read()
                complete_phrase_found = phrase in s
            file.close()
        if not complete_phrase_found:
            time.sleep(wait_time)

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
    
    run_info_batch = pd.read_csv(io.StringIO(results.stdout), sep='/t')
    run_info_batch = run_info_batch.loc[:,~run_info_batch.columns.str.startswith('Unnamed')]
    run_info_batch['JobID'] = run_info_batch['JobID'].str.replace('.batch' , '')
    run_info_batch['Max RAM Used (GB)'] = run_info_batch['MaxRSS'].str.replace('K' , '').astype(float)/1e6
    run_info_batch[['Allocated CPUs','Allocated RAM (GB)','Allocated Nodes']] = run_info_batch['AllocTRES'].str.split(',',expand=True)
    run_info_batch['Allocated CPUs'] = run_info_batch['Allocated CPUs'].str.replace('cpu=' , '').astype(int)
    run_info_batch['Allocated Nodes'] = run_info_batch['Allocated Nodes'].str.replace('node=' , '').astype(int)
    run_info_batch['Allocated RAM (GB)'] = run_info_batch['Allocated RAM (GB)'].str.replace('mem=' , '')
    run_info_batch['Allocated RAM (GB)'] = run_info_batch['Allocated RAM (GB)'].str.replace('G' , '').astype(float)
    run_info_batch['Elapsed'] = pd.to_timedelta(run_info_batch['Elapsed'].str.replace('-' , ' days '))
    run_info_batch['CPUTime'] = pd.to_timedelta(run_info_batch['CPUTime'].str.replace('-' , ' days '))
    colon_counts = run_info_batch['TotalCPU'].str.count(':')
    run_info_batch['TotalCPU'][colon_counts==1] = '00:' + run_info_batch['TotalCPU'][colon_counts==1]
    if any(run_info_batch['TotalCPU'].str.contains('-', regex=False)):
        run_info_batch['TotalCPU'].str.replace('-' , ' days ')
    run_info_batch['TotalCPU'] = pd.to_timedelta(run_info_batch['TotalCPU'])
    run_info_batch['CPU Efficiency (%)'] = 100*run_info_batch['TotalCPU']/run_info_batch['CPUTime']
    run_info_batch['RAM Efficiency (%)'] = 100*run_info_batch['Max RAM Used (GB)']/run_info_batch['Allocated RAM (GB)']
    
    job_ids_request = ','.join(job_ids)
    request = f"sacct --jobs={job_ids_request} --format=JobID,JobName,Timelimit -p --delimiter='/t'"
    results = run(request, shell=True, capture_output=True, text=True)
    run_info = pd.read_csv(io.StringIO(results.stdout), sep='/t')
    run_info = run_info[run_info['JobID'].isin(job_ids)]
    run_info = run_info.loc[:,~run_info.columns.str.startswith('Unnamed')]   
    run_info['Timelimit'] = pd.to_timedelta(run_info['Timelimit'].str.replace('-' , ' days '))
    
    
    job_stats = run_info.merge(run_info_batch, on='JobID')
    job_stats['Timelimit Used %'] =  100*job_stats['Elapsed']/job_stats['Timelimit']
    job_stats = job_stats[['JobID','JobName','Elapsed','Timelimit','Timelimit Used %', 'Allocated Nodes',
                         'Allocated CPUs', 'TotalCPU', 'CPUTime', 'CPU Efficiency (%)','Max RAM Used (GB)',
                         'Allocated RAM (GB)', 'RAM Efficiency (%)']]
    return job_stats