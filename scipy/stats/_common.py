
from dataclasses import dataclass


@dataclass
class ConfidenceInterval:
    """
    Class for confidence intervals.
    """
    low: float
    high: float
