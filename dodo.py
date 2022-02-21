"""
Info:
--------------------------
Command to list the tasks -
doit list

Command to run all the tasks -
doit

Checking all the info for a task-
doit info <task_name>

Targeted run for individual task -
doit doc_build
"""
# TODO: move the variables to dedicated config file
# TODO: improve doc messaging

from doit.task import clean_targets

DOIT_CONFIG = {'verbosity': 2}


def init_task():
    print("initiating build")


def enter_module_to_test():
    name_module = input("Enter the module to test: ")
    return name_module


def task_scipy_build():
    """
    Scipy build actions
    """
    name_module = enter_module_to_test()
    task_list = ["python dev.py --build-only", f"python dev.py --no-build -s % {name_module}"]
    for task_ in task_list:
        task_name = task_.split(' ')[2]
        yield{'name': task_name,
              'doc': 'initial build task',
              'actions': [task_],
              'file_dep': ["dev.py"],
              'clean': [clean_targets, init_task]
              }


def task_doc_build():
    """
    Runs document build tasks
    """
    return {'actions': ["python dev.py --doc"],
            'file_dep': ["dev.py"],
            'doc': 'doc build task'
            }
