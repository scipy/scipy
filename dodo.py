"""
Info: Run tests, builds and other tasks using doit
--------------------------
Command to list the tasks-
doit list

Command to check all the info for a task-
doit info <task_name>

Help command-
doit help <task name>

Targeted run for individual task-
Examples:
        $ doit build
        $ doit test -f <module name>
        $ doit doc-build
"""

from doit import create_after

DOIT_CONFIG = {'verbosity': 2}


# def paver_write_release_task():
#     """
#     Write Release task contents from pavement.py
#     """
#    # TODO: create a new script or add contents here


def task_build():
    """
    Scipy build task
    """
    return {'actions': ["python dev.py --build-only"],
            'file_dep': ["dev.py"],
            'doc': 'Task: Initializing build task'
            }


def task_test():
    """
    Runs the tests for a given module
    """
    return {'actions': ["python dev.py --no-build -s %(module)s"],
            'file_dep': ["dev.py"],
            'doc': 'Task: Initializing tests for the chosen module',
            'params': [{'name': 'module',
                        'short': 'f',
                        'default': '',
                        'type': str,
                        'help': 'Enter the module name to run tests'}],
            }


def task_doc_compile():
    """
    Runs document build tasks
    """
    return {'actions': ["python dev.py --doc"],
            'basename': 'doc-compile',
            'file_dep': ["dev.py"],
            'doc': 'Task: Initializing document build task'
            }


"""
Task dependency group
"""


def task_doc_build():
    """
    Task group with dependency for document build
    """
    return {'actions': None,
            'basename': 'doc-build',
            'doc': 'Task Group: Initializing document build tasks',
            'task_dep': ['build', 'doc-compile']}


def gen_release_tasks():
    """
    Task generator for release tasks
    """
    yield {'actions': ["python tools/authors.py"],
           'basename': 'release-authors',
           'file_dep': ["tools/authors.py"],
           'doc': 'Task: Initializing create author list'}
    yield {'actions': ["paver write_release_and_log"],
           'basename': 'release-notes',
           'file_dep': ["pavement.py"],
           'doc': 'Task: Initializing create release notes'}


def task_release():
    """
    Call to task generator
    """
    yield gen_release_tasks()
