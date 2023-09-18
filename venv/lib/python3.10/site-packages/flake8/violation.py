"""Contains the Violation error class used internally."""
import functools
import linecache
import logging
from typing import Dict
from typing import Match
from typing import NamedTuple
from typing import Optional
from typing import Set

from flake8 import defaults
from flake8 import utils


LOG = logging.getLogger(__name__)


@functools.lru_cache(maxsize=512)
def _find_noqa(physical_line: str) -> Optional[Match[str]]:
    return defaults.NOQA_INLINE_REGEXP.search(physical_line)


class Violation(NamedTuple):
    """Class representing a violation reported by Flake8."""

    code: str
    filename: str
    line_number: int
    column_number: int
    text: str
    physical_line: Optional[str]

    def is_inline_ignored(self, disable_noqa: bool) -> bool:
        """Determine if a comment has been added to ignore this line.

        :param disable_noqa:
            Whether or not users have provided ``--disable-noqa``.
        :returns:
            True if error is ignored in-line, False otherwise.
        """
        physical_line = self.physical_line
        # TODO(sigmavirus24): Determine how to handle stdin with linecache
        if disable_noqa:
            return False

        if physical_line is None:
            physical_line = linecache.getline(self.filename, self.line_number)
        noqa_match = _find_noqa(physical_line)
        if noqa_match is None:
            LOG.debug("%r is not inline ignored", self)
            return False

        codes_str = noqa_match.groupdict()["codes"]
        if codes_str is None:
            LOG.debug("%r is ignored by a blanket ``# noqa``", self)
            return True

        codes = set(utils.parse_comma_separated_list(codes_str))
        if self.code in codes or self.code.startswith(tuple(codes)):
            LOG.debug(
                "%r is ignored specifically inline with ``# noqa: %s``",
                self,
                codes_str,
            )
            return True

        LOG.debug(
            "%r is not ignored inline with ``# noqa: %s``", self, codes_str
        )
        return False

    def is_in(self, diff: Dict[str, Set[int]]) -> bool:
        """Determine if the violation is included in a diff's line ranges.

        This function relies on the parsed data added via
        :meth:`~StyleGuide.add_diff_ranges`. If that has not been called and
        we are not evaluating files in a diff, then this will always return
        True. If there are diff ranges, then this will return True if the
        line number in the error falls inside one of the ranges for the file
        (and assuming the file is part of the diff data). If there are diff
        ranges, this will return False if the file is not part of the diff
        data or the line number of the error is not in any of the ranges of
        the diff.

        :returns:
            True if there is no diff or if the error is in the diff's line
            number ranges. False if the error's line number falls outside
            the diff's line number ranges.
        """
        if not diff:
            return True

        # NOTE(sigmavirus24): The parsed diff will be a defaultdict with
        # a set as the default value (if we have received it from
        # flake8.utils.parse_unified_diff). In that case ranges below
        # could be an empty set (which is False-y) or if someone else
        # is using this API, it could be None. If we could guarantee one
        # or the other, we would check for it more explicitly.
        line_numbers = diff.get(self.filename)
        if not line_numbers:
            return False

        return self.line_number in line_numbers
