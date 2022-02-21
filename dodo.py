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
        $ doit doc
"""

from doit.task import clean_targets

DOIT_CONFIG = {'verbosity': 2}


def init_task():
    print("initiating scipy tasks")


def task_build():
    """
    Scipy build task
    """
    return {'actions': ["python dev.py --build-only"],
            'file_dep': ["dev.py"],
            'doc': 'Initializing build task'
            }


def task_test():
    """
    Runs the tests for a given module
    """
    return {'actions': ["python dev.py --no-build -s %(module)s"],
            'file_dep': ["dev.py"],
            'doc': 'Initializing tests for the chosen module',
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
            'file_dep': ["dev.py"],
            'doc': 'Initializing document build task'
            }


def task_doc():
    return {'actions': None,
            'doc': 'Initializing document build tasks',
            'task_dep': ['build', 'doc_compile']}
