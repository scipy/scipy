#!/usr/bin/env python3
"""
A dummy script that only echos its input arguments.

This is useful in case a platform-independent way to run a no-op command
on a target in a meson.build file is needed (e.g., to establish a
dependency between targets).
"""
import logging
import sys

logging.debug(f"Passed args to `scipy/_build_utils/echo.py`: {sys.argv[1:]}")
