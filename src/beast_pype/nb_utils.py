import papermill as pm
import nbformat as nbf


def make_kernelspec(kernel_name):
    """Create a kernelspec metadata dict for a Jupyter notebook.

    For Python kernels the name and display_name are the conda environment
    name.  For R kernels they are the conda environment name appended with
    ``_R``.  For bash kernels they are ``bash``.

    Parameters
    ----------
    kernel_name : str
        Name of the kernel.  If it ends with ``_R``, the language is set
        to ``'R'``.  If it equals ``'bash'``, the language is ``'bash'``.
        Otherwise the language is ``'python'``.

    Returns
    -------
    dict
        A kernelspec dict with ``'name'``, ``'display_name'`` and
        ``'language'`` keys.
    """
    if kernel_name.lower() == 'bash':
        language = 'bash'
    elif kernel_name.endswith('_R'):
        language = 'R'
    else:
        language = 'python'

    return {
        'name': kernel_name,
        'display_name': kernel_name,
        'language': language,
    }


def change_notebook_kernel(input_path, output_path, kernel):
    """
    Change a Jupyter Notebook's kernel_name.

    Parameters
    ----------
    input_path: str
        Path to the input notebook.
    output_path: str
        Path to the output notebook.
    kernel: str
        Name of the kernel_name to change to.
    """
    notebook = nbf.read(input_path, as_version=4)
    notebook['metadata']['kernelspec'] = make_kernelspec(kernel)
    nbf.write(notebook, output_path)

def execute_notebook(input_path,
                     output_path,
                     kernel_name=None,
                     **kwargs
                     ):
    """
    Execute a Jupyter Notebook using a specific Jupyter kernel_name.

    Notes
    ------
    This function is wrapper for papermill.execute, as at papermill v 2.6.0 the
    `kernal_name` argument for that function is not working.

    Parameters
    ----------
    input_path: str
        Path to input notebook.
    output_path: str
        Path to save executed notebook.
    kernel_name: str, optional
        Name of the kernel_name to execute.
    kwargs:
        Additional arguments to pass to papermill.execute.
    """
    if kernel_name is not None:
        change_notebook_kernel(input_path, output_path, kernel_name)
        pm.execute_notebook(input_path=output_path,
                            output_path=output_path,
                            kernel_name=kernel_name, **kwargs)
    else:
        pm.execute_notebook(input_path=input_path, output_path=output_path, **kwargs)