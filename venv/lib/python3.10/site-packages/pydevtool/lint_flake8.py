from flake8.main.application import Application as Flake8App
from flake8.checker import FileChecker
from flake8.formatting.default import Default as _Flake8Formatter
from doit.exceptions import TaskFailed


class Flake8Formatter(_Flake8Formatter):
    def write(self, line, source):
        self.console.print(line, source)


# FIXME: this is way slower compared to using flake8 CLI directly.
# flake8 uses all processors available, doit multiprocessing is failing due non-pickable
# thread.RLock reference. rich?
class LintFlake8():
    """flake8"""
    def __init__(self, console, config_file=None):
        self.console = console
        self.app = Flake8App()
        self.app.initialize([])
        self.checks = self.app.file_checker_manager.checks.to_dictionary()

        # It was supposed to be a plugin
        # https://flake8.pycqa.org/en/latest/plugin-development/index.html
        formatter = Flake8Formatter(self.app.options)
        formatter.console = console
        self.guide = self.app.file_checker_manager.style_guide.default_style_guide
        self.guide.formatter = formatter

    def excluded(self, fn):
        return self.app.file_checker_manager.is_path_excluded(fn)

    def __call__(self, fn):
        """execute pyflakes and pycodestyle on single file"""
        manager = self.app.file_checker_manager
        checker = FileChecker(str(fn), self.checks, manager.options)
        filename, results, stats = checker.run_checks()
        success = True
        for (error_code, line_number, column, text, physical_line) in results:
            has_reported_error = self.guide.handle_error(
                code=error_code,
                filename=filename,
                line_number=line_number,
                column_number=column,
                text=text,
                physical_line=physical_line,
            )
            if has_reported_error:
                success = False
        if(not success):
            return TaskFailed('Check failed', report=False)
