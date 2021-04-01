"""
Functions for analysis of variance
==================================

.. currentmodule:: scipy.stats.anova

.. autosummary::
   :toctree: generated/

   oneway
   oneway_from_labels_values
   twoway_from_data_grid
   twoway_from_a_b_values
   twoway_nestes_from_data_grid
   twoway_nested_from_a_b_values

"""

from ._anova_oneway import oneway, oneway_from_labels_values
from ._anova_twoway import (twoway_from_a_b_values,
                            twoway_from_data_grid)
from ._anova_twoway_nested import (twoway_nested_from_data_grid,
                                   twoway_nested_from_a_b_values)
