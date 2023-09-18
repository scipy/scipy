from pathlib import Path

import pycodestyle


class CodeStyleRichReporter(pycodestyle.BaseReport):
    console = None  # must be assigned after creation

    def __init__(self, options):
        super().__init__(options)
        self._fmt = pycodestyle.REPORT_FORMAT.get(
            options.format.lower(), options.format)

    def error(self, line_number, offset, text, check):
        """Report an error, according to options."""
        code = text[:4]
        if self._ignore_code(code):
            return
        if code in self.counters:
            self.counters[code] += 1
        else:
            self.counters[code] = 1
            self.messages[code] = text[5:]
        # Don't care about expected errors or warnings
        if code in self.expected:
            return
        if self.print_filename and not self.file_errors:
            print(self.filename)
        self.file_errors += 1
        self.total_errors += 1

        self.console.print(self._fmt % {
            'path': self.filename,
            'row': self.line_offset + line_number, 'col': offset + 1,
            'code': code, 'text': text[5:],
        })
        return code


class LintCodeStyle():
    """pycodestyle"""
    def __init__(self, console, config_file=None):
        style = pycodestyle.StyleGuide(config_file=config_file)
        style_reporter = style.init_report(CodeStyleRichReporter)
        style_reporter.console = console
        self.style = style
        # FIXME: hack to convert exlude dirs into patterns
        self.style.options.exclude = [
            p if p.endswith('.py') else f'{p}/*'
            for p in self.style.options.exclude]


    def iter_paths(self, paths):
        for path in paths:
            for match in Path(path).rglob('*.py'):
                fname = str(match)
                if not self.style.excluded(fname):
                    yield fname

    def __call__(self, fn):
        """execute pyflakes and pycodestyle on single file"""
        style_result = self.style.input_file(fn)
        return style_result == 0
