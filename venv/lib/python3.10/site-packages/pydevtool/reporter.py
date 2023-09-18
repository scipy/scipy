"""doit reporter using rich, supports progress bar"""

from rich.console import Console
from rich.progress import Progress, BarColumn, TimeElapsedColumn
from rich.theme import Theme
from rich.panel import Panel
from doit.exceptions import UnmetDependency, TaskError, TaskFailed as DoitTaskFailed


custom_theme = Theme({
    "name": "dim cyan",
    "event": "yellow",
    "failure": "red",
    "success": "green",
    "up-to-date": "cyan",
})


class DoitDebugReporter():
    desc = 'debug reporter - print out reporter events'

    def __init__(self, outstream, options):
        self.console = Console(theme=custom_theme)

    def initialize(self, tasks, selected_tasks):
        """called just after tasks have been loaded before execution starts"""
        self.console.print('[event]INITIALIZE [name]')
        # console.print(Panel(str(tasks), title="tasks"))
        # console.print(Panel(str(selected_tasks), title="selected"))

    def get_status(self, task):
        self.console.print(f'[event]STATUS [name]{task.name}')

    def execute_task(self, task):
        self.console.print(f'[event]EXECUTE [name]{task.name}')

    def skip_uptodate(self, task):
        self.console.print(f'[up-to-date]UP-TO-DATE [name]{task.name}')

    def skip_ignore(self, task):
        self.console.print(f'[event]SKIP [name]{task.name}')

    def add_failure(self, task, exception):
        self.console.print(f'[failure]FAILURE [name]{task.name}')
        if not isinstance(exception, DoitTaskFailed):
            self.console.print(Panel(str(exception), title="Error"))

    def add_success(self, task):
        self.console.print(f'[success]SUCCESS [name]{task.name}')

    def complete_run(self):
        self.console.print('[event]COMPLETE [name]')

    def cleanup_error(self, exception):
        self.console.print('[event]CLEANUP ERROR [name]')

    def runtime_error(self, msg):
        self.console.print('[event]ERROR [name]')
        self.console.print(Panel(msg, title="Error"))

    def teardown_task(self, task):
        """called when starts the execution of teardown action"""
        self.console.print('[event]TEARDOWN [name]')


class DoitRichReporter():
    desc = 'Rich color reporter'

    def __init__(self, console: Console, options=None):
        self.console = console
        columns = [
            "[progress.description]{task.description}",
            TimeElapsedColumn(),
            BarColumn(),
            "[progress.percentage]{task.completed} / {task.total}",
            "[progress.description]{task.fields[current]}",
        ]
        self.prog_bars = {}
        self.progress = Progress(*columns, console=self.console,
                                 redirect_stdout=False, redirect_stderr=False)


    ########## Progress bar management
    def add_progress_bar(self, description, total):
        # click API calls each progress bar a task
        # Add an extra field that contains the detail of the item (doit-task name)
        self.prog_bars[description] = self.progress.add_task(
            description, total=total, current=None)

    def update_progress(self, task, **kwargs):
        # set "current" on progress bar based on task name
        try:
            base, name = task.name.split(':', 1)
        except ValueError:
            # group task has not ':'
            if task.name in self.prog_bars:
                self.progress.update(self.prog_bars[task.name], current='')
            return
        if base in self.prog_bars:
            self.progress.update(self.prog_bars[base], current=name, **kwargs)

    ###############################


    def initialize(self, tasks, selected_tasks):
        """called just after tasks have been loaded before execution starts"""
        # console.print(Panel(str(tasks), title="tasks"))
        # console.print(Panel(str(selected_tasks), title="selected"))
        self.progress.start()

    def get_status(self, task):
        self.update_progress(task)

    def execute_task(self, task):
        pass

    def skip_uptodate(self, task):
        self.update_progress(task, advance=1)

    def skip_ignore(self, task):
        self.update_progress(task, advance=1)

    def add_failure(self, task, fail_info):
        if task.has_subtask and isinstance(fail_info, UnmetDependency):
            return
        if task.subtask_of:
            self.update_progress(task, advance=1)
        if fail_info.report:
            if fail_info.traceback:
                self.console.print(Panel(
                    "".join(fail_info.traceback),
                    title=f"{task.name}",
                    subtitle=fail_info.message,
                    border_style="red",
                ))
            else:
                label = 'Error' if isinstance(fail_info, TaskError) else 'Failed'
                self.console.print(f'[red]Task {label} - {task.name}'
                                   f' => {fail_info.message}')

    def add_success(self, task):
        self.update_progress(task, advance=1)

    def complete_run(self):
        self.progress.stop()

    def cleanup_error(self, exception):
        raise NotImplementedError()

    def runtime_error(self, msg):
        self.console.print(Panel(msg, title="Error"))
        # console.print("[red bold] msg")

    def teardown_task(self, task):
        """called when starts the execution of teardown action"""
        pass
