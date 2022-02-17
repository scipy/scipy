"""
Command to run the tasks -
doit -f <filename.py>

Targeted run -
doit -f updated_cli.py doc_build

Command to list the tasks -
doit -f updated_cli.py list
"""
from doit.task import clean_targets

"""
Setting config variable 
TODO: move the variables to dedicated config file 
"""
DOIT_CONFIG = {'verbosity': 2}


def init_task():
    print("initiating build")


def make_task(func):
    """make decorated function a task-creator"""
    func.create_doit_tasks = func
    return func


@make_task
def task_build():
    """
    Runs the build action with options (Hardcoded currently)
    """
    return {'actions': ["python dev.py --build-only"],
            'file_dep': ["dev.py"],
            'doc': 'initial build task',
            'clean': [clean_targets, init_task],
            }


@make_task
def task_test_action():
    """
    Runs cluster test tasks
    """
    return {'actions': ["python dev.py --no-build -s cluster"],
            'file_dep': ["dev.py"],
            'doc': 'cluster test task'
            }


@make_task
def task_doc_build():
    """
    Runs document build tasks
    """
    return {'actions': ["python dev.py --doc"],
            'file_dep': ["dev.py"],
            'doc': 'doc build task'
            }
