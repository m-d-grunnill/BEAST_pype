import importlib.resources as importlib_resources


def get_lib_resource_path(package, resource):
    """Get the path to a resource in a package.

    Parameters:
    package: str
        The name of the package containing the resource.
    resource: str
        The name of the resource.

    Returns:
    resource_path: pathlib.PosixPath
        Path to the resource.
    """
    resource_path, = importlib_resources.path(package, resource).args
    return resource_path

def path_to_workflows():
    """Get the path to the workflows directory.

    Returns:
    workflows_path: pathlib.PosixPath
        Path to the workflows directory.
    """
    workflows_path = get_lib_resource_path('beast_pype', 'workflows')
    return workflows_path

def path_to_report_templates():
    """Get the path to the report templates directory.

    Returns:
    reports_path: pathlib.PosixPath
        Path to the report templates directory.
    """
    reports_path = get_lib_resource_path('beast_pype', 'report_templates')
    return reports_path

def path_to_workflow_modules():
    """Get the path to the workflow modules directory.

    Returns:
    workflow_modules_path: pathlib.PosixPath
        Path to the workflow modules directory.
    """
    workflow_modules_path = get_lib_resource_path('beast_pype', 'workflow_modules')
    return workflow_modules_path